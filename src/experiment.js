// Experiments: draw many data sets from a model that is known, run a
// procedure on each, and compare what it says with the truth. Every result
// here is a simulated answer: it carries its seed, so it can be repeated, and
// it reports its own sampling error.
import { Rational, RationalInterval } from "@ratmath/core";
import { rng as makeRng } from "@ratmath/discrete";
import { refuse, toQ, Q, whole, probability, show, lowOf, highOf, isQ, toNumber } from "./num.js";
import { Model } from "./model.js";
import { TestResult, proportionInterval } from "./tests.js";

/** How each run gets its data: from { model, n }, or from draw(rng, i). */
function source({ model, n, draw }) {
  if (draw) {
    if (typeof draw !== "function") refuse("bad-argument", "draw must be a function of the generator");
    return draw;
  }
  if (!(model instanceof Model)) refuse("bad-argument", "an experiment needs { model, n } or a draw(rng) function");
  const size = whole(n, "the sample size n", { min: 1 });
  return (rng) => model.sample(size, rng);
}

/** A rate estimated by simulation: hits out of runs, with an exact interval for the rate it estimates. */
export class SimulatedRate {
  constructor(fields) {
    Object.assign(this, fields);
    this.kind = "simulated";
    Object.freeze(this);
  }
  /** Where a stated value lies against the interval: "below", "inside" or "above" it. */
  compare(value) {
    const q = toQ(value, "the value");
    return q.compareTo(this.interval.low) < 0 ? "below" : q.compareTo(this.interval.high) > 0 ? "above" : "inside";
  }
  toString() {
    return `${show(this.rate, 4)} (simulated: ${this.hits} of ${this.runs} runs, seed ${JSON.stringify(this.seed)}; ${show(this.interval.level.multiply(Q(100)), 1)}% interval ${show(this.interval.low, 4)} to ${show(this.interval.high, 4)})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/**
 * How often something happens: `event(data)` says yes or no for each of
 * `runs` data sets.
 *
 *   howOften({ model, n, event, runs: 1000, seed: 1 })
 *   howOften({ draw: (rng) => data, event, … })
 *
 * The result has `rate` = hits/runs exactly, and `interval`, the exact
 * (Clopper–Pearson) interval at `level` for the true rate. Run i uses the
 * stream rng(seed).fork("run i"), so runs do not depend on each other or on
 * how many there are.
 */
export function howOften({ event, runs = 1000, seed = 1, level = "0.95", ...rest }) {
  if (typeof event !== "function") refuse("bad-argument", "howOften needs an event: a function of the data that returns true or false");
  const draw = source(rest);
  const R = whole(runs, "runs", { min: 1 });
  const base = makeRng(seed);
  let hits = 0;
  for (let i = 0; i < R; i++) if (event(draw(base.fork(`run ${i}`), i), i)) hits++;
  const ci = proportionInterval(hits, R, { level });
  return new SimulatedRate({ rate: Q(hits, R), hits, runs: R, seed, interval: { low: ci.low, high: ci.high, level: ci.level } });
}

/**
 * How often an interval procedure contains the truth.
 *   coverage({ model, n, interval: (data) => [low, high], truth, runs, seed })
 * `interval` may return [low, high] or { low, high }, as numbers, Rationals
 * or enclosures. The result is a SimulatedRate with `meanWidth` added.
 */
export function coverage({ interval, truth, ...rest }) {
  if (typeof interval !== "function") refuse("bad-argument", "coverage needs an interval procedure: a function of the data that returns [low, high]");
  if (truth === undefined) refuse("bad-argument", "coverage needs the truth the intervals are meant to contain");
  const t = toQ(truth, "the truth");
  let width = 0;
  const result = howOften({
    ...rest,
    event(data) {
      const got = interval(data);
      const [lo, hi] = Array.isArray(got) ? got : [got.low, got.high];
      const [l, h] = [toQ(lowOf(asExact(lo)), "the lower end"), toQ(highOf(asExact(hi)), "the upper end")];
      width += h.subtract(l).toNumber();
      return l.compareTo(t) <= 0 && t.compareTo(h) <= 0;
    },
  });
  return new SimulatedRate({ ...result, truth: t, meanWidth: width / result.runs });
}
const asExact = (v) => (v instanceof Rational || v instanceof RationalInterval ? v : toQ(v, "an end of the interval"));

/** A p-value or a yes/no from whatever a test procedure returns. */
function rejects(outcome, alpha) {
  if (typeof outcome === "boolean") return outcome;
  const p = outcome instanceof TestResult ? outcome.pValue : outcome;
  if (typeof p === "number") return p <= alpha.toNumber();
  if (isQ(p) || p instanceof RationalInterval) {
    if (highOf(p).compareTo(alpha) <= 0) return true;
    if (lowOf(p).compareTo(alpha) > 0) return false;
    return toNumber(p) <= alpha.toNumber(); // an enclosure that straddles α: decided by its midpoint
  }
  return refuse("bad-argument", "a test procedure must return a p-value, a test result, or true/false");
}
/**
 * How often a test rejects at level α.
 *   rejectionRate({ model, n, test: (data) => pValue, alpha: 0.05, runs, seed })
 * When the model satisfies the null hypothesis this estimates the test's
 * true size; otherwise its power.
 */
export function rejectionRate({ test, alpha = "0.05", ...rest }) {
  if (typeof test !== "function") refuse("bad-argument", "rejectionRate needs a test: a function of the data that returns a p-value");
  const a = probability(alpha, "α");
  const result = howOften({ ...rest, event: (data) => rejects(test(data), a) });
  return new SimulatedRate({ ...result, alpha: a });
}

/**
 * The rejection rate at each of several effects.
 *   powerCurve({ effects: [0, 0.2, 0.5], model: (effect) => model(...), n, test, alpha, runs, seed })
 * or with draw: (effect) => (rng) => data. Every effect uses the same seed,
 * so the curve is smooth in the sense that the same random numbers are
 * reused. Returns [{ effect, rate, interval, … }].
 */
export function powerCurve({ effects, model, draw, ...rest }) {
  if (!Array.isArray(effects) || !effects.length) refuse("bad-argument", "powerCurve needs a list of effects");
  if (typeof model !== "function" && typeof draw !== "function") refuse("bad-argument", "powerCurve needs model: (effect) => model or draw: (effect) => (rng) => data");
  return effects.map((effect) => {
    const r = rejectionRate({ ...rest, ...(typeof model === "function" ? { model: model(effect) } : { draw: draw(effect) }) });
    return new SimulatedRate({ ...r, effect });
  });
}

/** What an estimator did over many data sets, in floating point. */
export class SimulatedEstimates {
  constructor(fields) {
    Object.assign(this, fields);
    this.kind = "simulated";
    Object.freeze(this);
  }
  toString() {
    const f = (v) => Number(v.toPrecision(5));
    const bias = this.truth === undefined ? "" : `, bias ${f(this.bias)} ± ${f(1.96 * this.standardError)}, rmse ${f(this.rmse)}`;
    return `mean ${f(this.mean)}, sd ${f(this.sd)}${bias} (simulated: ${this.runs} runs, seed ${JSON.stringify(this.seed)})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}
/**
 * The sampling distribution of an estimator.
 *   estimatorStudy({ model, n, estimator: (data) => number, truth, runs, seed })
 * Returns the estimates and their mean, sd and standardError (of the mean);
 * with a truth, also bias and rmse. bias ± 1.96·standardError is an
 * approximate 95% interval for the true bias.
 */
export function estimatorStudy({ estimator, truth, runs = 1000, seed = 1, ...rest }) {
  if (typeof estimator !== "function") refuse("bad-argument", "estimatorStudy needs an estimator: a function of the data that returns a number");
  const draw = source(rest);
  const R = whole(runs, "runs", { min: 2 });
  const base = makeRng(seed);
  const estimates = new Array(R);
  for (let i = 0; i < R; i++) {
    const v = estimator(draw(base.fork(`run ${i}`), i), i);
    estimates[i] = typeof v === "number" ? v : toNumber(v);
  }
  const mean = estimates.reduce((s, v) => s + v, 0) / R;
  const sd = Math.sqrt(estimates.reduce((s, v) => s + (v - mean) ** 2, 0) / (R - 1));
  const out = { estimates, mean, sd, standardError: sd / Math.sqrt(R), runs: R, seed };
  if (truth !== undefined) {
    const t = typeof truth === "number" ? truth : toNumber(toQ(truth, "the truth"));
    out.truth = t;
    out.bias = mean - t;
    out.rmse = Math.sqrt(estimates.reduce((s, v) => s + (v - t) ** 2, 0) / R);
  }
  return new SimulatedEstimates(out);
}
