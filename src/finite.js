// Finite probability, exactly: a distribution is a list of outcomes with
// rational probabilities that add up to 1. Events, conditioning, random
// variables, joint and marginal distributions and Bayes' rule are all sums
// and quotients of those rationals.
import { Rational, RationalInterval, Integer } from "@ratmath/core";
import { combinations, permutations, words } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, Q, sqrt, whole, probability } from "./num.js";

/**
 * A text key that is the same for outcomes that should count as the same and
 * different for outcomes that differ. Equal numbers agree whatever their type
 * (1, 1n, Integer 1 and Rational 1/1 all key as `q1/1`). Keys of tuple and set
 * members are length-prefixed, so a member's text can never run into the next
 * member's: ["a,sb", "c"] and ["a", "b,sc"] stay different outcomes.
 */
export function keyOf(o) {
  if (o instanceof Rational) return `q${o.numerator}/${o.denominator}`;
  if (o instanceof Integer) return `q${o.value}/1`;
  switch (typeof o) {
    case "number":
      return Number.isSafeInteger(o) ? `q${o}/1` : Number.isFinite(o) ? keyOf(toQ(o)) : `n${o}`;
    case "bigint":
      return `q${o}/1`;
    case "string":
      return `s${o}`;
    case "boolean":
      return `b${o}`;
    case "undefined":
      return "u";
    default:
      if (o === null) return "null";
      if (Array.isArray(o)) return `[${o.map(member).join("")}]`;
      if (o instanceof RationalInterval) return `i${keyOf(o.low)}:${keyOf(o.high)}`;
      if (o instanceof Set) return `{${[...new Set([...o].map(member))].sort().join("")}}`;
      return `o${JSON.stringify(o, objectValue)}`;
  }
}
// A member's key with its length in front: "3:sab". Joined members can then be
// read back one at a time, so different lists of members give different keys.
const member = (o) => {
  const k = keyOf(o);
  return `${k.length}:${k}`;
};
// Inside plain objects, keep values apart by type the way keyOf does at the
// top level: the string "1/2" and the rational 1/2, or "5n" and 5n, differ,
// while 5, 5n and Rational 5 agree.
// The holder's raw value is used because JSON.stringify applies toJSON first.
function objectValue(k, v) {
  const raw = k === "" ? v : this[k];
  return typeof raw === "string" || typeof raw === "bigint" || typeof raw === "number" || typeof raw === "boolean" || raw instanceof Rational || raw instanceof Integer || raw instanceof RationalInterval
    ? keyOf(raw)
    : v;
}
const numeric = (o, what) => {
  if (o instanceof Rational) return o;
  if (o instanceof Integer || typeof o === "number" || typeof o === "bigint") return toQ(o, what);
  if (typeof o === "boolean") return o ? ONE : ZERO;
  return refuse("not-numeric", `${what} needs numbers, but an outcome is ${showOutcome(o)}; use map(f) to turn outcomes into numbers`);
};
export const showOutcome = (o) => (Array.isArray(o) ? `(${o.map(showOutcome).join(", ")})` : o instanceof Set ? `{${[...o].map(showOutcome).join(", ")}}` : String(o));
/** The most outcomes a distribution will list at once; products of independent parts are kept implicit and may be far larger. */
export const LIMIT = 2000000;
const PRODUCT = Symbol("product");

// One factor of an implicit product, laid out for streaming: its outcomes in
// order, each with an integer weight over the common denominator `den`, and
// an index from outcome key to position.
function layout(dist) {
  const entries = dist.entries;
  const den = entries.reduce((l, [, p]) => (l * p.denominator) / gcd(l, p.denominator), 1n);
  const outs = entries.map((e) => e[0]);
  const w = entries.map(([, p]) => (p.numerator * den) / p.denominator);
  const index = new Map(outs.map((o, i) => [keyOf(o), i]));
  return Object.freeze({ dist, outs, w, den, index, uniform: w.every((x) => x === 1n) });
}

/**
 * A finite distribution: outcomes with exact rational probabilities that add
 * up to 1. It is held in one of two ways.
 *
 *   a table     every outcome and its probability, listed (finite, map, given, …)
 *   a product   independent parts, kept implicit (product, power, dice)
 *
 * A product of k parts has as many outcomes as the sizes of the parts
 * multiplied, so it is never listed unless asked: one outcome's probability is
 * the product of the parts' probabilities, an event's probability is counted
 * by streaming the outcomes one at a time with integer weights and dividing
 * once, map(f) keeps only the distinct values of f, marginals are the parts
 * themselves, and a sample is drawn part by part. Every answer is the same as
 * for the listed table, outcome order and seeded samples included.
 */
export class Finite {
  #index;
  #parts = null; // the parts of a product, laid out
  #count = 0n; // the number of outcomes of a product
  #den = 1n; // the common denominator of a product's outcome probabilities
  #uniform = false; // every outcome of the product equally likely
  #listed = null; // a product's entries, once asked for
  /** entries: [outcome, probability] pairs; equal outcomes are merged. The probabilities must add up to 1. */
  constructor(entries, { normalize = false, [PRODUCT]: parts } = {}) {
    if (parts) {
      this.#parts = Object.freeze(parts.map(layout));
      this.#count = this.#parts.reduce((n, f) => n * BigInt(f.outs.length), 1n);
      this.#den = this.#parts.reduce((n, f) => n * f.den, 1n);
      this.#uniform = this.#parts.every((f) => f.uniform);
      Object.freeze(this);
      return;
    }
    const index = new Map();
    const list = [];
    for (const [o, w] of entries) {
      const p = w instanceof Rational ? w : toQ(w, "a probability");
      if (p.numerator < 0n) refuse("bad-argument", `a probability cannot be negative (${p} for ${showOutcome(o)})`);
      if (p.numerator === 0n) continue;
      const k = keyOf(o);
      const at = index.get(k);
      if (at === undefined) {
        index.set(k, list.length);
        list.push([o, p]);
      } else list[at][1] = list[at][1].add(p);
    }
    let total = ZERO;
    for (const [, p] of list) total = total.add(p);
    if (total.numerator === 0n) refuse("bad-argument", "a distribution needs at least one outcome with positive probability");
    if (!total.equals(ONE)) {
      if (!normalize) refuse("bad-argument", `the probabilities add up to ${total}, not 1; pass { normalize: true } to rescale them`);
      for (const e of list) e[1] = e[1].divide(total);
    }
    Object.defineProperty(this, "entries", { value: Object.freeze(list.map((e) => Object.freeze(e))), enumerable: true });
    this.#index = index;
    Object.freeze(this);
  }
  /**
   * Every [outcome, probability] pair, in order. A table has them already; a
   * product lists them on first use, and refuses beyond LIMIT outcomes.
   */
  get entries() {
    if (!this.#listed) {
      this.#tooMany("list every outcome");
      const list = [];
      this.#scan((t, w) => void list.push(Object.freeze([t, this.#q(w)])));
      this.#listed = Object.freeze(list);
    }
    return this.#listed;
  }
  /** The number of outcomes, as a number. */
  get size() {
    return this.#parts ? Number(this.#count) : this.entries.length;
  }
  /** The number of outcomes, exactly. */
  get count() {
    return this.#parts ? this.#count : BigInt(this.entries.length);
  }
  /** For a product, the independent parts (one Finite per coordinate); null for a table. */
  get factors() {
    return this.#parts ? this.#parts.map((f) => f.dist) : null;
  }
  get outcomes() {
    return this.entries.map((e) => e[0]);
  }
  /** Each [outcome, probability] pair in turn, without listing them all; works for products of any size. */
  *each() {
    if (!this.#parts) return yield* this.entries;
    const k = this.#parts.length;
    const d = new Array(k).fill(0);
    for (;;) {
      let w = 1n;
      const t = this.#parts.map((f, i) => ((w *= f.w[d[i]]), f.outs[d[i]]));
      yield [t, this.#q(w)];
      let i = k - 1;
      while (i >= 0 && ++d[i] === this.#parts[i].outs.length) d[i--] = 0;
      if (i < 0) return;
    }
  }
  #q(w) {
    return new Rational(w, this.#den);
  }
  #tooMany(what) {
    if (this.#count > BigInt(LIMIT))
      refuse(
        "too-large",
        `this product has ${this.#count} outcomes, too many to ${what} (the limit is ${LIMIT}); ask for what you need instead: p(outcome), prob(event), map(f), marginal(i), sumOf(n), sample(rng) and each() all work without listing the outcomes`,
      );
  }
  /** A numeric question asked of tuples: refuse it the way a table would, without listing anything. */
  #notNumeric(what) {
    return numeric(this.#parts.map((f) => f.outs[0]), what);
  }
  /**
   * Visit every outcome of a product in order: visit(tuple, weight), where the
   * weight is an integer over #den (1n for every outcome when all are equally
   * likely). The tuple is a fresh array each time. Returning true stops.
   */
  #scan(visit) {
    const parts = this.#parts;
    const k = parts.length;
    const d = new Array(k).fill(0);
    const t = parts.map((f) => f.outs[0]);
    const uniform = this.#uniform;
    const pw = new Array(k);
    const prefix = (from) => {
      for (let j = from; j < k; j++) pw[j] = (j ? pw[j - 1] : 1n) * parts[j].w[d[j]];
    };
    if (!uniform) prefix(0);
    for (;;) {
      if (visit(t.slice(), uniform ? 1n : pw[k - 1])) return;
      let i = k - 1;
      while (i >= 0 && ++d[i] === parts[i].outs.length) {
        d[i] = 0;
        t[i] = parts[i].outs[0];
        i--;
      }
      if (i < 0) return;
      t[i] = parts[i].outs[d[i]];
      if (!uniform) prefix(i);
    }
  }
  /** The total integer weight of the outcomes where holds(tuple) is true. */
  #weigh(holds) {
    if (this.#uniform) {
      let n = 0;
      this.#scan((t) => {
        if (holds(t)) n++;
      });
      return BigInt(n);
    }
    let s = 0n;
    this.#scan((t, w) => {
      if (holds(t)) s += w;
    });
    return s;
  }
  /** Where a tuple sits in the product's order, or -1 when it is not an outcome. */
  #position(o) {
    if (!Array.isArray(o) || o.length !== this.#parts.length) return -1n;
    let r = 0n;
    for (let i = 0; i < o.length; i++) {
      const f = this.#parts[i];
      const at = f.index.get(keyOf(o[i]));
      if (at === undefined) return -1n;
      r = r * BigInt(f.outs.length) + BigInt(at);
    }
    return r;
  }
  #holds(event) {
    if (typeof event === "function") return event;
    if (event instanceof Set) {
      const keys = new Set([...event].map(keyOf));
      return (o) => keys.has(keyOf(o));
    }
    const k = keyOf(event);
    return (o) => keyOf(o) === k;
  }
  /** The distinct outcomes of a Set or single-outcome event, in the product's order. */
  #members(event) {
    const seen = new Map();
    for (const o of event instanceof Set ? event : [event]) {
      const at = this.#position(o);
      if (at >= 0n && !seen.has(at)) seen.set(at, this.#parts.map((f, i) => f.outs[f.index.get(keyOf(o[i]))]));
    }
    return [...seen.entries()].sort((a, b) => (a[0] < b[0] ? -1 : 1)).map((e) => e[1]);
  }
  /** The probability of one outcome. */
  p(outcome) {
    if (this.#parts) {
      if (this.#position(outcome) < 0n) return ZERO;
      return this.#parts.reduce((s, f, i) => s.multiply(f.dist.p(outcome[i])), ONE);
    }
    const at = this.#index.get(keyOf(outcome));
    return at === undefined ? ZERO : this.entries[at][1];
  }
  /**
   * The probability of an event: a function of the outcome that says yes or
   * no, a Set of outcomes, or one outcome. With a second event, the
   * conditional probability P(event | given).
   */
  prob(event, given) {
    if (given !== undefined && this.#parts) {
      // Count both events in one pass instead of listing the condition.
      const g = this.#holds(given);
      const e = this.#holds(event);
      let both = 0n;
      let total = 0n;
      this.#scan((t, w) => {
        if (!g(t)) return;
        total += w;
        if (e(t)) both += w;
      });
      if (total === 0n) refuse("impossible-event", "cannot condition on an event of probability 0");
      return new Rational(both, total);
    }
    if (given !== undefined) return this.given(given).prob(event);
    if (typeof event !== "function" && !(event instanceof Set)) return this.p(event);
    if (this.#parts) {
      if (event instanceof Set) return this.#members(event).reduce((s, o) => s.add(this.p(o)), ZERO);
      return this.#q(this.#weigh(event));
    }
    const holds = this.#holds(event);
    let s = ZERO;
    for (const [o, p] of this.entries) if (holds(o)) s = s.add(p);
    return s;
  }
  /** The distribution conditioned on an event of positive probability. */
  given(event) {
    if (this.#parts) {
      let kept;
      if (typeof event === "function") {
        kept = [];
        this.#scan((t, w) => {
          if (!event(t)) return;
          if (kept.length === LIMIT) refuse("too-large", `more than ${LIMIT} outcomes satisfy the condition; ask prob(event, given) for a conditional probability instead`);
          kept.push([t, this.#q(w)]);
        });
      } else kept = this.#members(event).map((o) => [o, this.p(o)]);
      if (!kept.length) refuse("impossible-event", "cannot condition on an event of probability 0");
      return new Finite(kept, { normalize: true });
    }
    const holds = this.#holds(event);
    const kept = this.entries.filter(([o]) => holds(o));
    if (!kept.length) refuse("impossible-event", "cannot condition on an event of probability 0");
    return new Finite(kept, { normalize: true });
  }
  /** The distribution of f(outcome): a random variable. A product is streamed, keeping only the distinct values. */
  map(f) {
    if (!this.#parts) return new Finite(this.entries.map(([o, p]) => [f(o), p]));
    const index = new Map();
    const rows = [];
    this.#scan((t, w) => {
      const v = f(t);
      const k = keyOf(v);
      const at = index.get(k);
      if (at !== undefined) rows[at][1] += w;
      else {
        if (rows.length === LIMIT) refuse("too-large", `the function takes more than ${LIMIT} different values`);
        index.set(k, rows.length);
        rows.push([v, w]);
      }
    });
    return new Finite(rows.map(([v, w]) => [v, this.#q(w)]));
  }
  /** Keep only some outcomes and rescale: the same as given(event). */
  filter(event) {
    return this.given(event);
  }
  /**
   * The joint distribution of independent draws from this and the others:
   * outcomes are tuples. Nothing is listed; the result keeps the parts.
   */
  product(...others) {
    return new Finite(null, { [PRODUCT]: [this, ...others.map(asFinite)] });
  }
  /** n independent copies, as tuples. */
  power(n) {
    const k = whole(n, "the number of copies", { min: 1 });
    return this.product(...Array(k - 1).fill(this));
  }
  /** The distribution of X + Y for an independent Y (numbers only). */
  plus(other) {
    if (this.#parts) this.#notNumeric("plus");
    const d = asFinite(other);
    const out = [];
    for (const [a, p] of this.entries) for (const [b, q] of d.entries) out.push([numeric(a, "plus").add(numeric(b, "plus")), p.multiply(q)]);
    return new Finite(out);
  }
  /** The distribution of the sum of n independent copies (numbers only), by repeated doubling. */
  sumOf(n) {
    if (this.#parts) this.#notNumeric("sumOf");
    let k = whole(n, "the number of copies", { min: 1 });
    let base = this.map((o) => numeric(o, "sumOf"));
    let acc = null;
    while (k > 0) {
      if (k & 1) acc = acc ? acc.plus(base) : base;
      k >>= 1;
      if (k) base = base.plus(base);
    }
    return acc;
  }
  /**
   * For tuple outcomes: the distribution of the sum of the coordinates. A
   * product of numeric parts convolves the parts (no outcome is visited), so
   * dice(20).sum() is immediate; a table of tuples adds up each outcome.
   */
  sum() {
    if (this.#parts) return this.#parts.map((f) => f.dist.map((o) => numeric(o, "sum"))).reduce((a, b) => a.plus(b));
    return this.map((o) => {
      if (!Array.isArray(o)) refuse("bad-argument", "sum needs outcomes that are tuples; for the sum of n independent copies use sumOf(n)");
      return o.reduce((s, x) => s.add(numeric(x, "sum")), ZERO);
    });
  }
  /** For tuple outcomes: the distribution of the chosen coordinates. For a product these are its parts. */
  marginal(...which) {
    const parts = this.#parts;
    if (parts && which.length && which.every((i, j) => Number.isInteger(i) && i >= 0 && i < parts.length && (j === 0 || i > which[j - 1]))) {
      if (which.length === 1) return parts[which[0]].dist;
      return new Finite(null, { [PRODUCT]: which.map((i) => parts[i].dist) });
    }
    return this.map((o) => {
      if (!Array.isArray(o)) refuse("bad-argument", "marginal needs outcomes that are tuples");
      return which.length === 1 ? o[which[0]] : which.map((i) => o[i]);
    });
  }
  /** The joint distribution of several functions of the outcome. */
  joint(...fs) {
    return this.map((o) => fs.map((f) => f(o)));
  }
  /**
   * Are two events independent, P(A and B) = P(A)·P(B)? Given two functions
   * that return values rather than yes or no, tests whether the two random
   * variables are independent: every pair of values must factor.
   */
  independent(a, b) {
    const fa = typeof a === "function" ? a : this.#holds(a);
    const fb = typeof b === "function" ? b : this.#holds(b);
    const j = this.joint(fa, fb);
    const ma = j.marginal(0);
    const mb = j.marginal(1);
    for (const [x, px] of ma.entries) for (const [y, py] of mb.entries) if (!j.p([x, y]).equals(px.multiply(py))) return false;
    return true;
  }
  /** E[f(X)], or E[X] for numeric outcomes. */
  mean(f) {
    if (this.#parts) return f ? this.map((o) => numeric(f(o), "mean")).mean() : this.#notNumeric("mean");
    let s = ZERO;
    for (const [o, p] of this.entries) s = s.add(numeric(f ? f(o) : o, "mean").multiply(p));
    return s;
  }
  expect(f) {
    return this.mean(f);
  }
  /** The k-th moment E[Xᵏ]; with { central: true }, E[(X − μ)ᵏ]. */
  moment(k, { central = false, of } = {}) {
    if (this.#parts) return of ? this.map(of).moment(k, { central }) : this.#notNumeric("moment");
    const n = whole(k, "the order of the moment");
    const mu = central ? this.mean(of) : ZERO;
    return this.mean((o) => numeric(of ? of(o) : o, "moment").subtract(mu).pow(n));
  }
  variance(f) {
    return this.moment(2, { central: true, of: f });
  }
  /** The standard deviation: exact when the variance is a perfect square, otherwise a certified enclosure. */
  sd(f, { bits = 96 } = {}) {
    return sqrt(this.variance(f), bits);
  }
  covariance(f, g) {
    if (this.#parts) return this.joint(f, g).covariance((o) => o[0], (o) => o[1]);
    const mf = this.mean(f);
    const mg = this.mean(g);
    return this.mean((o) => numeric(f(o), "covariance").subtract(mf).multiply(numeric(g(o), "covariance").subtract(mg)));
  }
  /** The correlation of f and g: exact when it is rational, otherwise certified. */
  correlation(f, g, { bits = 96 } = {}) {
    if (this.#parts) return this.joint(f, g).correlation((o) => o[0], (o) => o[1], { bits });
    const c = this.covariance(f, g);
    const vv = this.variance(f).multiply(this.variance(g));
    if (vv.numerator === 0n) refuse("undefined", "the correlation is not defined when one variable is constant");
    const r = sqrt(c.multiply(c).divide(vv), bits);
    if (c.numerator >= 0n) return r;
    return r instanceof Rational ? r.negate() : new RationalInterval(r.high.negate(), r.low.negate());
  }
  /** The numeric outcomes in increasing order with their probabilities. */
  get support() {
    if (this.#parts) this.#notNumeric("support");
    return this.entries.map(([o, p]) => [numeric(o, "support"), p]).sort((a, b) => a[0].compareTo(b[0]));
  }
  /** P(X ≤ x). */
  cdf(x) {
    const q = toQ(x, "x");
    if (this.#parts) this.#notNumeric("cdf");
    let s = ZERO;
    for (const [o, p] of this.entries) if (numeric(o, "cdf").compareTo(q) <= 0) s = s.add(p);
    return s;
  }
  /** The smallest value x with P(X ≤ x) ≥ p. */
  quantile(p) {
    const level = probability(p, "the level");
    let s = ZERO;
    const sup = this.support;
    for (const [x, px] of sup) {
      s = s.add(px);
      if (s.compareTo(level) >= 0) return x;
    }
    /* c8 ignore next */
    return sup[sup.length - 1][0];
  }
  median() {
    return this.quantile(Q(1, 2));
  }
  /** The outcomes of greatest probability. For a product, every combination of the parts' modes. */
  modes() {
    if (this.#parts) {
      const each = this.#parts.map((f) => f.dist.modes());
      const n = each.reduce((m, ms) => m * BigInt(ms.length), 1n);
      if (n > BigInt(LIMIT)) refuse("too-large", `there are ${n} most likely outcomes, more than the limit of ${LIMIT}; ask each part with factors[i].modes()`);
      return each.reduce((acc, ms) => acc.flatMap((t) => ms.map((m) => [...t, m])), [[]]);
    }
    let best = ZERO;
    for (const [, p] of this.entries) if (p.compareTo(best) > 0) best = p;
    return this.entries.filter(([, p]) => p.equals(best)).map((e) => e[0]);
  }
  /**
   * One outcome drawn with the exact probabilities, from a seeded generator.
   * A product draws one integer below the common denominator and reads the
   * parts off it in turn, the same outcome its listed table would give.
   */
  sample(rng) {
    if (this.#parts) {
      let r = rng.int(this.#den);
      let rest = this.#den;
      return this.#parts.map((f) => {
        rest /= f.den;
        const s = r / rest;
        let cum = 0n;
        let j = 0;
        while (s >= cum + f.w[j]) cum += f.w[j++];
        r = ((s - cum) * rest + (r % rest)) / f.w[j];
        return f.outs[j];
      });
    }
    const den = this.entries.reduce((l, [, p]) => (l * p.denominator) / gcd(l, p.denominator), 1n);
    let r = rng.int(den);
    for (const [o, p] of this.entries) {
      const w = (p.numerator * den) / p.denominator;
      if (r < w) return o;
      r -= w;
    }
    /* c8 ignore next */
    return this.entries[this.entries.length - 1][0];
  }
  samples(n, rng) {
    return Array.from({ length: n }, () => this.sample(rng));
  }
  /** Rows { outcome, probability }, in the order given or sorted by a comparison. */
  table({ sort } = {}) {
    const rows = this.entries.map(([outcome, probability]) => ({ outcome, probability }));
    if (sort === true) rows.sort((a, b) => numeric(a.outcome, "sort").compareTo(numeric(b.outcome, "sort")));
    else if (typeof sort === "function") rows.sort((a, b) => sort(a.outcome, b.outcome));
    return rows;
  }
  toString() {
    const first = [];
    if (this.#parts) this.#scan((t, w) => first.push([t, this.#q(w)]) === 12);
    else first.push(...this.entries.slice(0, 12));
    const shown = first.map(([o, p]) => `${showOutcome(o)}: ${p}`);
    return `{ ${shown.join(", ")}${this.count > 12n ? `, … (${this.count} outcomes)` : ""} }`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return `Finite ${this.toString()}`;
  }
}
const gcd = (a, b) => {
  while (b) [a, b] = [b, a % b];
  return a;
};
const isFamily = (v) => v && typeof v === "object" && typeof v.unrank === "function" && typeof v.count === "bigint";
const asFinite = (v) => (v instanceof Finite ? v : finite(v));

/**
 * A finite distribution from
 *   an array of outcomes          equally likely; a repeated outcome counts each time
 *   an array and weights          finite(["a", "b"], ["1/3", "2/3"])
 *   an object or a Map            finite({ H: "1/2", T: "1/2" })
 *   a family from @ratmath/discrete   every object of the family equally likely
 * Weights that do not add up to 1 are refused unless { normalize: true }.
 */
export function finite(spec, weights, options) {
  if (spec instanceof Finite) return spec;
  if (weights && !Array.isArray(weights)) [weights, options] = [undefined, weights];
  if (isFamily(spec)) {
    if (spec.count > BigInt(LIMIT)) refuse("too-large", `the family has ${spec.count} objects; the limit is ${LIMIT}`);
    const n = Number(spec.count);
    const p = Q(1, n);
    return new Finite(Array.from(spec, (o) => [o, p]));
  }
  if (spec instanceof Map) return new Finite([...spec.entries()], options);
  if (Array.isArray(spec)) {
    if (!spec.length) refuse("bad-argument", "a distribution needs at least one outcome");
    if (weights) {
      if (weights.length !== spec.length) refuse("bad-argument", "one weight for each outcome");
      return new Finite(spec.map((o, i) => [o, weights[i]]), options);
    }
    const p = Q(1, spec.length);
    return new Finite(spec.map((o) => [o, p]));
  }
  if (spec && typeof spec === "object") return new Finite(Object.entries(spec), options);
  return refuse("bad-argument", "finite needs an array of outcomes, an object or Map of probabilities, or a family");
}

/** A fair die: the numbers 1 to sides. */
export const die = (sides = 6) => finite(Array.from({ length: whole(sides, "the number of sides", { min: 1 }) }, (_, i) => i + 1));
/** n dice thrown together: outcomes are tuples. */
export const dice = (n, sides = 6) => die(sides).power(n);
/** A coin showing "H" with probability p. */
export const coin = (p = "1/2") => {
  const q = probability(p, "the probability of heads");
  return new Finite([["H", q], ["T", ONE.subtract(q)]]);
};
/** A standard deck of 52 cards, as text such as "Q♥". */
export const deck = () => {
  const cards = [];
  for (const suit of ["♠", "♥", "♦", "♣"]) for (const rank of ["A", "2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K"]) cards.push(rank + suit);
  return cards;
};
/**
 * Drawing k of the items, every result equally likely.
 *   { replace: false, ordered: false }   a hand: k distinct items, order ignored (the default)
 *   { ordered: true }                    k distinct items in the order drawn
 *   { replace: true, ordered: true }     k draws with the item put back each time
 *   { replace: true }                    not offered: unordered draws with replacement are not equally likely
 */
export function draw(items, k, { replace = false, ordered = false } = {}) {
  const n = whole(k, "the number drawn");
  if (replace && !ordered) refuse("bad-argument", "unordered draws with replacement are not equally likely; use { ordered: true } and map the result");
  return finite(replace ? words(items, n) : ordered ? permutations(items, n) : combinations(items, n));
}

/**
 * Bayes' rule on a finite set of hypotheses. `prior` is a distribution over
 * hypotheses; `likelihood(h)` is the probability of what was observed if h is
 * true (or an object or Map from hypothesis to that probability).
 *
 *   { posterior, evidence, table: [{ hypothesis, prior, likelihood, joint, posterior }] }
 */
export function bayes(prior, likelihood) {
  const d = asFinite(prior);
  const like =
    typeof likelihood === "function"
      ? likelihood
      : likelihood instanceof Map
        ? (h) => [...likelihood.entries()].find(([k]) => keyOf(k) === keyOf(h))?.[1]
        : (h) => likelihood[h];
  const rows = d.entries.map(([hypothesis, p]) => {
    const raw = like(hypothesis);
    if (raw === undefined) refuse("bad-argument", `no likelihood was given for the hypothesis ${showOutcome(hypothesis)}`);
    const l = probability(raw, "a likelihood");
    return { hypothesis, prior: p, likelihood: l, joint: p.multiply(l) };
  });
  const evidence = rows.reduce((s, r) => s.add(r.joint), ZERO);
  if (evidence.numerator === 0n) refuse("impossible-event", "what was observed has probability 0 under every hypothesis");
  for (const r of rows) r.posterior = r.joint.divide(evidence);
  return { posterior: new Finite(rows.map((r) => [r.hypothesis, r.posterior])), evidence, table: rows };
}

