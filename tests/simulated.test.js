import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import * as S from "../index.js";

describe("models", () => {
  const line = S.model({ truth: [2, 3], x: S.uniform(0, 10), noise: S.normal(0, 1) });
  test("a sample remembers where it came from and repeats from its seed", () => {
    const d = line.sample(50, 42);
    expect(d.kind).toBe("simulated");
    expect(d.n).toBe(50);
    expect(d.seed).toBe(42);
    expect(d.x.length).toBe(50);
    expect(d.y).toEqual(line.sample(50, 42).y);
    expect(d.y).not.toEqual(line.sample(50, 43).y);
    expect(d.trueY[0]).toBeCloseTo(2 + 3 * d.trueX[0], 12);
    expect(String(d)).toContain("50 simulated observations");
    expect(String(line)).toBe("model(truth 2 + 3·x; x from uniform(0, 10); noise normal(0, 1))");
    const fit = S.fast.leastSquares(d.x, d.y);
    expect(Math.abs(fit.slope - 3)).toBeLessThan(0.2);
    // a generator may be passed in place of a seed
    expect(line.sample(5, S.rng(42)).y).toEqual(line.sample(5, 42).y);
  });
  test("changing one part leaves the other parts' numbers alone", () => {
    const a = line.sample(30, 7);
    const b = line.with({ noise: S.normal(0, 5) }).sample(30, 7);
    expect(b.x).toEqual(a.x);
    const c = line.with({ contamination: { rate: "1/5", noise: S.cauchy(0, 10) } }).sample(30, 7);
    expect(c.x).toEqual(a.x);
    const dirty = c.contaminated.filter(Boolean).length;
    expect(dirty).toBeGreaterThan(0);
    expect(dirty).toBeLessThan(15);
    // the clean observations draw their noise from the same stream
    expect(c.y[c.contaminated.indexOf(false)]).toBeDefined();
  });
  test("every departure", () => {
    const m = S.model({
      truth: (x) => 1 + x * x, x: [0, 1, 2, 3], noise: S.normal(0, 1), spread: (x) => 0.1 * (1 + x),
      xNoise: S.normal(0, "0.05"), round: { x: "0.1", y: "0.5" }, censor: { below: 2, above: 8 }, missing: "1/4",
      contamination: { rate: "1/4", value: S.uniform(100, 200) },
    });
    const d = m.sample(40, 3);
    expect(d.trueX.slice(0, 5)).toEqual([0, 1, 2, 3, 0]);
    expect(d.trueY.slice(0, 4)).toEqual([1, 2, 5, 10]);
    expect(d.censored.some(Boolean)).toBe(true);
    expect(d.missing.some(Boolean)).toBe(true);
    expect(d.y.filter((v) => !Number.isNaN(v)).every((v) => v >= 2 && v <= 8 && Math.abs(v * 2 - Math.round(v * 2)) < 1e-9)).toBe(true);
    expect(d.x.every((v) => Math.abs(v * 10 - Math.round(v * 10)) < 1e-9)).toBe(true);
    const kept = d.complete();
    expect(kept.y.length).toBe(d.missing.filter((f) => !f).length);
    expect(kept.x.length).toBe(kept.y.length);
    expect(String(m)).toContain("censored below 2 and above 8");
    expect(String(m)).toContain("missing with probability 1/4");
    expect(String(m)).toContain("x recorded with error");
    expect(String(S.model({ truth: 5, name: "five" }))).toBe("five");
    expect(String(S.model({ truth: 5, round: "0.1", censor: { above: 9 }, noise: (rng) => rng.float() }))).toContain("rounded to 0.1");
  });
  test("a model with no predictor, and sources of every kind", () => {
    const m = S.model({ truth: 10, noise: S.finite([-1, 0, 1]) });
    const d = m.sample(200, 1);
    expect(d.x).toBeUndefined();
    expect(new Set(d.y)).toEqual(new Set([9, 10, 11]));
    expect(d.complete().x).toBeUndefined();
    expect(S.model({ truth: 0, noise: S.poisson(3) }).sample(5, 1).y.every(Number.isInteger)).toBe(true);
    expect(S.model({ truth: 0, noise: (rng) => rng.between(1, 6) }).sample(5, 1).y.every((v) => v >= 1 && v <= 6)).toBe(true);
    expect(S.model({ truth: 1, noise: "0.5" }).sample(2, 1).y).toEqual([1.5, 1.5]);
    expect(S.model({ truth: 7 }).sample(2, 1).y).toEqual([7, 7]);
    expect(S.model({ truth: [1, 2], x: (rng) => rng.between(0, 3) }).sample(6, 1).y.every((v) => [1, 3, 5, 7].includes(v))).toBe(true);
  });
  test("exact mode gives rationals that repeat everywhere", () => {
    const m = S.model({ truth: [2, 3], x: S.uniform(0, 10), noise: S.normal(0, 1), round: "0.01", censor: { below: 0, above: 40 }, missing: "1/10" });
    const d = m.sample(12, 42, { exact: true });
    expect(d.exact).toBe(true);
    expect(d.x.every((v) => v instanceof Rational)).toBe(true);
    expect(d.y.every((v) => v === null || v instanceof Rational)).toBe(true);
    expect(d.y.map(String)).toEqual(m.sample(12, 42, { exact: true }).y.map(String));
    const kept = d.complete();
    const fit = S.leastSquares(kept.x, kept.y);
    expect(fit.slope instanceof Rational).toBe(true);
    expect(Math.abs(fit.slope.toNumber() - 3)).toBeLessThan(0.6);
    const e = S.model({ truth: (x) => x.multiply(x), x: ["1/2", "1/3"], noise: S.finite([0]), spread: () => "2", xNoise: S.finite([1]) }).sample(2, 1, { exact: true });
    expect(e.y.map(String)).toEqual(["1/4", "1/9"]);
    expect(e.x.map(String)).toEqual(["3/2", "4/3"]);
    expect(S.model({ truth: "1/3" }).sample(1, 1, { exact: true }).y.map(String)).toEqual(["1/3"]);
    expect(S.model({ truth: 0, noise: S.binomial(3, "1/2") }).sample(3, 1, { exact: true }).y.every((v) => v instanceof Rational)).toBe(true);
    expect(S.model({ truth: 0, noise: (rng) => rng.between(1, 2) }).sample(3, 1, { exact: true }).y.every((v) => v instanceof Rational)).toBe(true);
  });
  test("mistakes are refused", () => {
    expect(() => S.model({ truth: 1, nose: S.normal() })).toThrow(/no part called "nose"/);
    expect(() => S.model({ noise: S.normal() })).toThrow(/needs a truth/);
    expect(() => S.model({ truth: 1, contamination: { rate: "1/2" } })).toThrow(/contamination needs/);
    expect(() => S.model({ truth: 1, contamination: { rate: 2, noise: S.normal() } })).toThrow(/between 0 and 1/);
    expect(() => S.model({ truth: 1, missing: -1 })).toThrow(/between 0 and 1/);
    expect(() => S.model({ truth: 1 }).sample(0)).toThrow(/at least 1/);
  });
});

describe("experiments", () => {
  test("coverage of intervals for a proportion", () => {
    const run = (method) =>
      S.coverage({
        draw: (rng) => S.binomial(30, "1/10").sample(rng),
        interval: (k) => {
          const ci = S.proportionInterval(k, 30, { method });
          return [ci.low, ci.high];
        },
        truth: "1/10", runs: 400, seed: 5,
      });
    const wald = run("wald");
    expect(wald.kind).toBe("simulated");
    expect(wald.seed).toBe(5);
    expect(wald.runs).toBe(400);
    expect(String(wald.rate)).toBe(String(new Rational(BigInt(wald.hits), 400n)));
    // the textbook interval covers far less often than it claims
    expect(wald.compare("0.95")).toBe("above");
    expect(wald.compare("0.5")).toBe("below");
    expect(wald.compare(wald.rate)).toBe("inside");
    const exact = run("exact");
    expect(exact.rate.toNumber()).toBeGreaterThanOrEqual(0.95);
    expect(exact.meanWidth).toBeGreaterThan(wald.meanWidth);
    expect(String(wald)).toMatch(/^0\.\d{4} \(simulated: \d+ of 400 runs, seed 5; 95% interval 0\.\d{4} to 0\.\d{4}\)$/);
    // the same seed gives the same answer
    expect(run("wald").hits).toBe(wald.hits);
    const obj = S.coverage({ draw: () => 0, interval: () => ({ low: -1, high: 1 }), truth: 0, runs: 10 });
    expect(String(obj.rate)).toBe("1");
    expect(() => S.coverage({ draw: () => 0, truth: 0 })).toThrow(/interval procedure/);
    expect(() => S.coverage({ draw: () => 0, interval: () => [0, 1] })).toThrow(/needs the truth/);
  });
  test("the size and power of a test", () => {
    const draw = (p) => (rng) => S.binomial(40, p).sample(rng);
    const size = S.rejectionRate({ draw: draw("1/2"), test: (k) => S.binomialTest(k, 40, "1/2"), runs: 500, seed: 11 });
    // an exact test never rejects more often than its level
    expect(size.interval.low.toNumber()).toBeLessThan(0.05);
    expect(String(size.alpha)).toBe("1/20");
    const curve = S.powerCurve({ effects: ["1/2", "13/20", "4/5"], draw, test: (k) => S.binomialTest(k, 40, "1/2").pValue, runs: 300, seed: 11 });
    expect(curve.map((c) => c.effect)).toEqual(["1/2", "13/20", "4/5"]);
    expect(curve[0].rate.toNumber()).toBeLessThan(curve[1].rate.toNumber());
    expect(curve[1].rate.toNumber()).toBeLessThan(curve[2].rate.toNumber());
    expect(curve[2].rate.toNumber()).toBeGreaterThan(0.9);
    const viaModel = S.powerCurve({ effects: [0, 2], model: (shift) => S.model({ truth: shift, noise: S.normal(0, 1) }), n: 10, test: (d) => S.fast.mean(d.y) > 1, runs: 100, seed: 2 });
    expect(viaModel[0].rate.toNumber()).toBeLessThan(0.1);
    expect(viaModel[1].rate.toNumber()).toBeGreaterThan(0.9);
    // p-values of every kind
    expect(String(S.rejectionRate({ draw: () => 0, test: () => 0.01, runs: 5 }).rate)).toBe("1");
    expect(String(S.rejectionRate({ draw: () => 0, test: () => S.poisson(1).pmf(9), runs: 5 }).rate)).toBe("1");
    expect(String(S.rejectionRate({ draw: () => 0, test: () => S.poisson(1).pmf(1), runs: 5 }).rate)).toBe("0");
    expect(() => S.rejectionRate({ draw: () => 0, test: () => "yes", runs: 5 })).toThrow(/must return a p-value/);
    expect(() => S.rejectionRate({ draw: () => 0, runs: 5 })).toThrow(/needs a test/);
    expect(() => S.powerCurve({ effects: [], draw, test: () => true })).toThrow(/list of effects/);
    expect(() => S.powerCurve({ effects: [1], test: () => true })).toThrow(/needs model/);
  });
  test("the bias of an estimator", () => {
    const m = S.model({ truth: 0, noise: S.normal(0, 2) });
    const divideByN = S.estimatorStudy({ model: m, n: 5, estimator: (d) => S.fast.variance(d.y, { population: true }), truth: 4, runs: 4000, seed: 1 });
    const divideByN1 = S.estimatorStudy({ model: m, n: 5, estimator: (d) => S.fast.variance(d.y), truth: 4, runs: 4000, seed: 1 });
    expect(divideByN.kind).toBe("simulated");
    // dividing by n is too small by the factor (n − 1)/n: the bias is −4/5
    expect(Math.abs(divideByN.bias + 0.8)).toBeLessThan(4 * divideByN.standardError);
    expect(Math.abs(divideByN1.bias)).toBeLessThan(4 * divideByN1.standardError);
    expect(divideByN.estimates.length).toBe(4000);
    expect(divideByN.rmse).toBeGreaterThan(0);
    expect(String(divideByN)).toMatch(/^mean [\d.]+, sd [\d.]+, bias -[\d.]+ ± [\d.]+, rmse [\d.]+ \(simulated: 4000 runs, seed 1\)$/);
    const noTruth = S.estimatorStudy({ model: m, n: 5, estimator: (d) => S.mean(d.y.map((v) => v.toFixed(2))), runs: 20 });
    expect(noTruth.bias).toBeUndefined();
    expect(String(noTruth)).not.toContain("bias");
    expect(() => S.estimatorStudy({ model: m, n: 5 })).toThrow(/needs an estimator/);
  });
  test("how often, in general", () => {
    const r = S.howOften({ model: S.model({ truth: 0, noise: S.normal(0, 1) }), n: 4, event: (d) => S.fast.max(d.y) > 2, runs: 600, seed: "max" });
    // P(max of 4 > 2) = 1 − Φ(2)⁴ = 0.0879…
    expect(r.compare("0.0879")).toBe("inside");
    expect(r.seed).toBe("max");
    expect(String(S.howOften({ draw: (rng, i) => i, event: (i) => i % 2 === 0, runs: 10, level: "0.9" }).rate)).toBe("1/2");
    expect(() => S.howOften({ draw: () => 0, runs: 5 })).toThrow(/needs an event/);
    expect(() => S.howOften({ event: () => true })).toThrow(/needs \{ model, n \}/);
    expect(() => S.howOften({ draw: 5, event: () => true })).toThrow(/function of the generator/);
    expect(() => S.howOften({ model: S.model({ truth: 0 }), event: () => true })).toThrow(/sample size/);
  });
});
