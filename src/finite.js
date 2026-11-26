// Finite probability, exactly: a distribution is a list of outcomes with
// rational probabilities that add up to 1. Events, conditioning, random
// variables, joint and marginal distributions and Bayes' rule are all sums
// and quotients of those rationals.
import { Rational, RationalInterval } from "@ratmath/core";
import { combinations, permutations, words } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, Q, sqrt, whole, probability } from "./num.js";

/** A text key that is the same for outcomes that should count as the same. */
export function keyOf(o) {
  if (o instanceof Rational) return `q${o.numerator}/${o.denominator}`;
  switch (typeof o) {
    case "number":
      return Number.isFinite(o) ? keyOf(toQ(o)) : `n${o}`;
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
      if (Array.isArray(o)) return `[${o.map(keyOf).join(",")}]`;
      if (o instanceof RationalInterval) return `i${o.low}:${o.high}`;
      if (o instanceof Set) return `{${[...o].map(keyOf).sort().join(",")}}`;
      return `o${JSON.stringify(o, (_k, v) => (typeof v === "bigint" ? `${v}n` : v instanceof Rational ? v.toString() : v))}`;
  }
}
const numeric = (o, what) => {
  if (o instanceof Rational) return o;
  if (typeof o === "number" || typeof o === "bigint") return toQ(o, what);
  if (typeof o === "boolean") return o ? ONE : ZERO;
  return refuse("not-numeric", `${what} needs numbers, but an outcome is ${showOutcome(o)}; use map(f) to turn outcomes into numbers`);
};
export const showOutcome = (o) => (Array.isArray(o) ? `(${o.map(showOutcome).join(", ")})` : o instanceof Set ? `{${[...o].map(showOutcome).join(", ")}}` : String(o));
const LIMIT = 2000000;

export class Finite {
  #index;
  /** entries: [outcome, probability] pairs; equal outcomes are merged. The probabilities must add up to 1. */
  constructor(entries, { normalize = false } = {}) {
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
    this.entries = Object.freeze(list.map((e) => Object.freeze(e)));
    this.#index = index;
    Object.freeze(this);
  }
  get size() {
    return this.entries.length;
  }
  get outcomes() {
    return this.entries.map((e) => e[0]);
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
  /** The probability of one outcome. */
  p(outcome) {
    const at = this.#index.get(keyOf(outcome));
    return at === undefined ? ZERO : this.entries[at][1];
  }
  /**
   * The probability of an event: a function of the outcome that says yes or
   * no, a Set of outcomes, or one outcome. With a second event, the
   * conditional probability P(event | given).
   */
  prob(event, given) {
    if (given !== undefined) return this.given(given).prob(event);
    if (typeof event !== "function" && !(event instanceof Set)) return this.p(event);
    const holds = this.#holds(event);
    let s = ZERO;
    for (const [o, p] of this.entries) if (holds(o)) s = s.add(p);
    return s;
  }
  /** The distribution conditioned on an event of positive probability. */
  given(event) {
    const holds = this.#holds(event);
    const kept = this.entries.filter(([o]) => holds(o));
    if (!kept.length) refuse("impossible-event", "cannot condition on an event of probability 0");
    return new Finite(kept, { normalize: true });
  }
  /** The distribution of f(outcome): a random variable. */
  map(f) {
    return new Finite(this.entries.map(([o, p]) => [f(o), p]));
  }
  /** Keep only some outcomes and rescale: the same as given(event). */
  filter(event) {
    return this.given(event);
  }
  /** The joint distribution of independent draws from this and the others: outcomes are tuples. */
  product(...others) {
    let acc = this.entries.map(([o, p]) => [[o], p]);
    for (const other of others) {
      const d = asFinite(other);
      if (acc.length * d.size > LIMIT) refuse("too-large", `the product would have ${acc.length * d.size} outcomes; the limit is ${LIMIT}`);
      const next = [];
      for (const [t, p] of acc) for (const [o, q] of d.entries) next.push([[...t, o], p.multiply(q)]);
      acc = next;
    }
    return new Finite(acc);
  }
  /** n independent copies, as tuples. */
  power(n) {
    const k = whole(n, "the number of copies", { min: 1 });
    return this.product(...Array(k - 1).fill(this));
  }
  /** The distribution of X + Y for an independent Y (numbers only). */
  plus(other) {
    const d = asFinite(other);
    const out = [];
    for (const [a, p] of this.entries) for (const [b, q] of d.entries) out.push([numeric(a, "plus").add(numeric(b, "plus")), p.multiply(q)]);
    return new Finite(out);
  }
  /** The distribution of the sum of n independent copies (numbers only), by repeated doubling. */
  sumOf(n) {
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
  /** For tuple outcomes: the distribution of the chosen coordinates. */
  marginal(...which) {
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
    let s = ZERO;
    for (const [o, p] of this.entries) s = s.add(numeric(f ? f(o) : o, "mean").multiply(p));
    return s;
  }
  expect(f) {
    return this.mean(f);
  }
  /** The k-th moment E[Xᵏ]; with { central: true }, E[(X − μ)ᵏ]. */
  moment(k, { central = false, of } = {}) {
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
    const mf = this.mean(f);
    const mg = this.mean(g);
    return this.mean((o) => numeric(f(o), "covariance").subtract(mf).multiply(numeric(g(o), "covariance").subtract(mg)));
  }
  /** The correlation of f and g: exact when it is rational, otherwise certified. */
  correlation(f, g, { bits = 96 } = {}) {
    const c = this.covariance(f, g);
    const vv = this.variance(f).multiply(this.variance(g));
    if (vv.numerator === 0n) refuse("undefined", "the correlation is not defined when one variable is constant");
    const r = sqrt(c.multiply(c).divide(vv), bits);
    if (c.numerator >= 0n) return r;
    return r instanceof Rational ? r.negate() : new RationalInterval(r.high.negate(), r.low.negate());
  }
  /** The numeric outcomes in increasing order with their probabilities. */
  get support() {
    return this.entries.map(([o, p]) => [numeric(o, "support"), p]).sort((a, b) => a[0].compareTo(b[0]));
  }
  /** P(X ≤ x). */
  cdf(x) {
    const q = toQ(x, "x");
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
  /** The outcomes of greatest probability. */
  modes() {
    let best = ZERO;
    for (const [, p] of this.entries) if (p.compareTo(best) > 0) best = p;
    return this.entries.filter(([, p]) => p.equals(best)).map((e) => e[0]);
  }
  /** One outcome drawn with the exact probabilities, from a seeded generator. */
  sample(rng) {
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
    const shown = this.entries.slice(0, 12).map(([o, p]) => `${showOutcome(o)}: ${p}`);
    return `{ ${shown.join(", ")}${this.size > 12 ? `, … (${this.size} outcomes)` : ""} }`;
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

