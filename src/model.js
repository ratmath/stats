// Data from a model that is known. A model is a truth plus a description of
// how observations depart from it; a sample from it is simulated data that
// remembers its seed and what was done to each observation, so that a
// procedure can be judged against the truth that produced the data.
import { Rational } from "@ratmath/core";
import { rng as makeRng, Random } from "@ratmath/discrete";
import { refuse, toQ, ZERO, whole, probability } from "./num.js";
import { Continuous } from "./continuous.js";
import { Discrete } from "./discrete-dists.js";
import { Finite } from "./finite.js";

const isDist = (d) => d instanceof Continuous || d instanceof Discrete || d instanceof Finite;
const asRng = (seed) => (seed instanceof Random ? seed : makeRng(seed));
const num = (v) => (typeof v === "number" ? v : v instanceof Rational ? v.toNumber() : Number(v));

/** One draw from a distribution, a function of the generator, or a constant. */
function drawFrom(source, rng, exact, bits) {
  if (source instanceof Continuous) return exact ? source.sampleExact(rng, { bits }) : source.sample(rng);
  if (source instanceof Discrete || source instanceof Finite) {
    const v = source.sample(rng);
    return exact ? toQ(v, "a draw") : num(v);
  }
  if (typeof source === "function") {
    const v = source(rng);
    return exact ? toQ(v, "a draw") : num(v);
  }
  return exact ? toQ(source, "a value") : num(toQ(source, "a value"));
}
const describeSource = (s) => (isDist(s) ? String(s) : typeof s === "function" ? "a function" : String(s));

/**
 * Simulated data. `y` are the observations; `x` the recorded predictor
 * values when the model has one. `kind` is always "simulated", `seed` and
 * `model` say where the data came from, and the rest records what happened:
 *   trueX, trueY     the values before any departure (the truth)
 *   contaminated, censored, missing    one flag for each observation
 */
export class Dataset {
  constructor(fields) {
    Object.assign(this, fields);
    this.kind = "simulated";
    Object.freeze(this);
  }
  get n() {
    return this.y.length;
  }
  /** The observations that are present, as { x, y } with matching positions. */
  complete() {
    const keep = this.missing.map((m) => !m);
    return { x: this.x ? this.x.filter((_, i) => keep[i]) : undefined, y: this.y.filter((_, i) => keep[i]) };
  }
  toString() {
    const tally = (flags) => flags.reduce((s, f) => s + (f ? 1 : 0), 0);
    return `${this.n} simulated observations from ${this.model} (seed ${JSON.stringify(this.seed)}; ${tally(this.contaminated)} contaminated, ${tally(this.censored)} censored, ${tally(this.missing)} missing)`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/**
 * A known truth and how observations depart from it.
 *
 *   truth          a number (the value every observation would have), the
 *                  coefficients [b0, b1, …] of b0 + b1·x + …, or a function of x
 *   x              where the predictor values come from: a list (reused in
 *                  order), a distribution, or a function of the generator.
 *                  Leave it out for a model with no predictor.
 *   noise          a distribution added to each observation
 *   spread         a function of x that multiplies the noise (unequal spread)
 *   contamination  { rate, noise }: with probability `rate` the noise for an
 *                  observation comes from this second distribution instead;
 *                  or { rate, value }: the whole observation is replaced
 *   xNoise         a distribution added to the recorded x (y uses the true x)
 *   round          a step: observations are rounded to its multiples
 *                  (or { x, y } steps)
 *   censor         { below, above }: values beyond a limit are recorded as the limit
 *   missing        the probability that an observation is lost (recorded as NaN, or null in exact mode)
 *
 * sample(n, seed) draws n observations in floating point. With { exact: true }
 * every draw is an exact rational (see sampleExact), so that exact statistics
 * of the data are exactly reproducible; that is much slower.
 */
export class Model {
  constructor(spec) {
    const known = new Set(["truth", "x", "noise", "spread", "contamination", "xNoise", "round", "censor", "missing", "name"]);
    for (const k of Object.keys(spec)) if (!known.has(k)) refuse("bad-argument", `a model has no part called "${k}"`);
    if (spec.truth === undefined) refuse("bad-argument", "a model needs a truth: a number, coefficients, or a function of x");
    if (spec.contamination) {
      const c = spec.contamination;
      probability(c.rate, "the contamination rate");
      if (c.noise === undefined && c.value === undefined) refuse("bad-argument", "contamination needs { rate, noise } or { rate, value }");
    }
    if (spec.missing !== undefined) probability(spec.missing, "the probability of a missing value");
    this.spec = Object.freeze({ ...spec });
    Object.freeze(this);
  }
  /** The same model with some parts changed: model.with({ noise: normal(0, 3) }). */
  with(changes) {
    return new Model({ ...this.spec, ...changes });
  }
  /** The truth at x: what an observation would be with no departure at all. */
  truthAt(x, exact = false) {
    const t = this.spec.truth;
    if (typeof t === "function") return exact ? toQ(t(x), "the truth") : num(t(x));
    if (Array.isArray(t)) {
      if (exact) {
        const q = toQ(x ?? 0, "x");
        return t.reduceRight((acc, c) => acc.multiply(q).add(toQ(c, "a coefficient")), ZERO);
      }
      return t.reduceRight((acc, c) => acc * x + num(toQ(c, "a coefficient")), 0);
    }
    return exact ? toQ(t, "the truth") : num(toQ(t, "the truth"));
  }
  sample(n, seed = 1, { exact = false, bits = 24 } = {}) {
    const count = whole(n, "the number of observations", { min: 1 });
    const s = this.spec;
    const base = asRng(seed);
    // one stream for each part, so that changing one part leaves the others' numbers alone
    const g = { x: base.fork("x"), noise: base.fork("noise"), contamination: base.fork("contamination"), xNoise: base.fork("x noise"), missing: base.fork("missing") };
    const zero = exact ? ZERO : 0;
    const plus = (a, b) => (exact ? a.add(b) : a + b);
    const times = (a, b) => (exact ? a.multiply(toQ(b, "the spread")) : a * num(b));
    const roundTo = (v, step) => {
      if (step === undefined) return v;
      if (exact) {
        const st = toQ(step, "the rounding step");
        return new Rational(v.divide(st).round()).multiply(st);
      }
      const st = num(toQ(step, "the rounding step"));
      return Math.round(v / st) * st;
    };
    const stepX = s.round !== null && typeof s.round === "object" && !(s.round instanceof Rational) ? s.round.x : undefined;
    const stepY = s.round !== null && typeof s.round === "object" && !(s.round instanceof Rational) ? s.round.y : s.round;
    const below = s.censor?.below !== undefined ? (exact ? toQ(s.censor.below, "the lower limit") : num(toQ(s.censor.below, "the lower limit"))) : undefined;
    const above = s.censor?.above !== undefined ? (exact ? toQ(s.censor.above, "the upper limit") : num(toQ(s.censor.above, "the upper limit"))) : undefined;
    const less = (a, b) => (exact ? a.compareTo(b) < 0 : a < b);
    const rate = s.contamination ? probability(s.contamination.rate, "the contamination rate") : null;
    const lost = s.missing !== undefined ? probability(s.missing, "the probability of a missing value") : null;

    const hasX = s.x !== undefined;
    const out = { x: hasX ? [] : undefined, y: [], trueX: hasX ? [] : undefined, trueY: [], contaminated: [], censored: [], missing: [] };
    for (let i = 0; i < count; i++) {
      let x;
      if (hasX) {
        x = Array.isArray(s.x) ? (exact ? toQ(s.x[i % s.x.length], "x") : num(toQ(s.x[i % s.x.length], "x"))) : drawFrom(s.x, g.x, exact, bits);
        out.trueX.push(x);
      }
      const truth = this.truthAt(x, exact);
      out.trueY.push(truth);
      const dirty = rate !== null && g.contamination.bool(rate);
      let y;
      if (dirty && s.contamination.value !== undefined) y = drawFrom(s.contamination.value, g.contamination, exact, bits);
      else {
        const source = dirty ? s.contamination.noise : s.noise;
        let e = source === undefined ? zero : drawFrom(source, g.noise, exact, bits);
        if (s.spread) e = times(e, s.spread(exact ? x : x));
        y = plus(truth, e);
      }
      y = roundTo(y, stepY);
      let cut = false;
      if (below !== undefined && less(y, below)) [y, cut] = [below, true];
      if (above !== undefined && less(above, y)) [y, cut] = [above, true];
      const gone = lost !== null && g.missing.bool(lost);
      out.y.push(gone ? (exact ? null : NaN) : y);
      out.contaminated.push(dirty);
      out.censored.push(cut);
      out.missing.push(gone);
      if (hasX) {
        let seen = s.xNoise === undefined ? x : plus(x, drawFrom(s.xNoise, g.xNoise, exact, bits));
        seen = roundTo(seen, stepX);
        out.x.push(seen);
      }
    }
    return new Dataset({ ...out, seed: seed instanceof Random ? seed.state : seed, exact, model: this.toString() });
  }
  toString() {
    const s = this.spec;
    if (s.name) return s.name;
    const t = typeof s.truth === "function" ? "f(x)" : Array.isArray(s.truth) ? s.truth.map((c, i) => (i === 0 ? `${c}` : i === 1 ? `${c}·x` : `${c}·x^${i}`)).join(" + ") : String(s.truth);
    const parts = [`truth ${t}`];
    if (s.x !== undefined) parts.push(`x from ${Array.isArray(s.x) ? `a list of ${s.x.length}` : describeSource(s.x)}`);
    if (s.noise !== undefined) parts.push(`noise ${describeSource(s.noise)}${s.spread ? " scaled by a function of x" : ""}`);
    if (s.contamination) parts.push(`contamination ${s.contamination.rate} from ${describeSource(s.contamination.noise ?? s.contamination.value)}`);
    if (s.xNoise !== undefined) parts.push(`x recorded with error ${describeSource(s.xNoise)}`);
    if (s.round !== undefined) parts.push(`rounded to ${typeof s.round === "object" && !(s.round instanceof Rational) ? JSON.stringify(s.round) : s.round}`);
    if (s.censor) parts.push(`censored ${s.censor.below !== undefined ? `below ${s.censor.below}` : ""}${s.censor.below !== undefined && s.censor.above !== undefined ? " and " : ""}${s.censor.above !== undefined ? `above ${s.censor.above}` : ""}`);
    if (s.missing !== undefined) parts.push(`missing with probability ${s.missing}`);
    return `model(${parts.join("; ")})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}
export const model = (spec) => new Model(spec);
