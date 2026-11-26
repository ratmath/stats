import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import * as S from "../index.js";

const Q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));
const lo = (v) => (v instanceof Rational ? v : v.low);
const hi = (v) => (v instanceof Rational ? v : v.high);
const mid = (v) => (lo(v).toNumber() + hi(v).toNumber()) / 2;
const width = (v) => hi(v).subtract(lo(v)).toNumber();
/** n draws from a distribution, rounded to hundredths: small exact rationals. */
const sample = (dist, n, rng) => Array.from({ length: n }, () => Q(Math.round(dist.sample(rng) * 100), 100));
const fmean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const fvar = (a) => {
  const m = fmean(a);
  return a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1);
};

const x = [5.1, 4.9, 5.6, 5.8, 6.0, 5.3, 4.7, 5.5];
const y = [4.4, 4.9, 5.0, 4.6, 5.2, 4.8, 4.5];

describe("t and z tests, against floating point and against each other", () => {
  test("one-sample t: exact t², certified p agreeing with the float distribution", () => {
    const r = S.tTest(x, { mu: 5 });
    const xf = x;
    const t = (fmean(xf) - 5) / Math.sqrt(fvar(xf) / xf.length);
    expect(r.statisticSquared.toNumber()).toBeCloseTo(t * t, 12);
    expect(mid(r.statistic)).toBeCloseTo(t, 12);
    expect(String(r.df)).toBe("7");
    const T = S.studentT(7);
    expect(mid(r.pValue)).toBeCloseTo(2 * (1 - T.fast.cdf(Math.abs(t))), 9);
    expect(width(r.pValue)).toBeLessThan(1e-15);
    expect(r.kind).toBe("certified");
    // the two one-sided p-values add up to 1
    const g = S.tTest(x, { mu: 5, alternative: "greater" });
    const l = S.tTest(x, { mu: 5, alternative: "less" });
    expect(mid(g.pValue) + mid(l.pValue)).toBeCloseTo(1, 14);
    expect(mid(g.pValue)).toBeCloseTo(mid(r.pValue) / 2, 14);
  });
  test("the t interval is the set of μ the test keeps", () => {
    const ci = S.tTest(x).interval();
    expect(S.tTest(x, { mu: lo(ci.low).subtract(Q(1, 1e6)) }).reject("0.05")).toBe(true);
    expect(S.tTest(x, { mu: ci.low.add(Q(1, 1e6)) }).reject("0.05")).toBe(false);
    expect(S.tTest(x, { mu: ci.high.subtract(Q(1, 1e6)) }).reject("0.05")).toBe(false);
    expect(S.tTest(x, { mu: ci.high.add(Q(1, 1e6)) }).reject("0.05")).toBe(true);
    const m = S.meanInterval(x);
    expect(m.low.equals(ci.low) && m.high.equals(ci.high)).toBe(true);
    const one = S.meanInterval(x, { alternative: "greater" });
    expect(one.high).toBe(null);
    expect(one.low.compareTo(ci.low)).toBeGreaterThan(0);
  });
  test("Welch's degrees of freedom and the pooled test", () => {
    const w = S.tTest(x, y);
    const [v1, v2] = [fvar(x) / x.length, fvar(y) / y.length];
    const df = (v1 + v2) ** 2 / (v1 ** 2 / (x.length - 1) + v2 ** 2 / (y.length - 1));
    expect(w.df.toNumber()).toBeCloseTo(df, 12);
    expect(w.statisticSquared.toNumber()).toBeCloseTo((fmean(x) - fmean(y)) ** 2 / (v1 + v2), 12);
    // a pooled two-sample t² is exactly the one-way ANOVA F for two groups
    const p = S.tTest(x, y, { equalVariance: true });
    const a = S.anova([x, y]);
    expect(p.statisticSquared.equals(a.statistic)).toBe(true);
    expect(mid(p.pValue)).toBeCloseTo(mid(a.pValue), 14);
  });
  test("the paired test is the one-sample test on the differences", () => {
    const xs = x.slice(0, 7);
    const p = S.tTest(xs, y, { paired: true });
    const d = S.tTest(xs.map((v, i) => S.toQ(v).subtract(S.toQ(y[i]))));
    expect(p.statisticSquared.equals(d.statisticSquared)).toBe(true);
    expect(p.pairs).toBe(7);
  });
  test("z test against Φ", () => {
    const r = S.zTest(x, { sigma: "1/2", mu: 5 });
    const z = (fmean(x) - 5) / (0.5 / Math.sqrt(8));
    expect(mid(r.statistic)).toBeCloseTo(z, 12);
    expect(mid(r.pValue)).toBeCloseTo(2 * (1 - S.normal(0, 1).fast.cdf(z)), 10);
    expect(S.zTest([1, 2, 3], { sigma: 1, mu: 2 }).pValue.equals(Q(1))).toBe(true);
    expect(() => S.zTest(x)).toThrow(/known standard deviation/);
  });
  test("refusals", () => {
    expect(() => S.tTest([1, 1, 1])).toThrow(/standard error is 0/);
    expect(() => S.tTest([1])).toThrow(/at least two/);
    expect(() => S.tTest(x, y, { paired: true })).toThrow(/same length/);
    expect(() => S.tTest(x, { alternative: "bigger" })).toThrow(/alternative/);
    expect(() => S.tTest([[1, 2], 3])).toThrow(/exact numbers/);
  });
});

describe("variances, proportions, chi-square, correlation", () => {
  test("variance interval and test agree, and the F test", () => {
    const ci = S.varianceInterval(x);
    const v = ci.estimate;
    expect(String(v)).toBe(String(S.variance(x)));
    expect(S.varianceTest(x, ci.low.subtract(Q(1, 1e9))).reject()).toBe(true);
    expect(S.varianceTest(x, ci.low.add(Q(1, 1e9))).reject()).toBe(false);
    expect(S.varianceTest(x, ci.high.subtract(Q(1, 1e9))).reject()).toBe(false);
    expect(S.varianceTest(x, ci.high.add(Q(1, 1e9))).reject()).toBe(true);
    expect(ci.sd.low.multiply(ci.sd.low).compareTo(ci.low)).toBeLessThanOrEqual(0);
    const f = S.fTest(x, y);
    expect(f.statistic.equals(S.variance(x).divide(S.variance(y)))).toBe(true);
    const fi = f.interval();
    expect(S.fTest(x, y, { ratio: fi.low.subtract(Q(1, 1e9)) }).reject()).toBe(true);
    expect(S.fTest(x, y, { ratio: fi.high.subtract(Q(1, 1e9)) }).reject()).toBe(false);
  });
  test("proportions: one sample, two samples, and the 2×2 identity X² = z²", () => {
    const r = S.proportionTest(60, 100);
    expect(String(r.statisticSquared)).toBe("4");
    expect(mid(r.pValue)).toBeCloseTo(0.0455002638963584, 13);
    const c = S.proportionTest(60, 100, "1/2", { correct: true });
    expect(String(c.statisticSquared)).toBe("361/100");
    const two = S.proportionTest([60, 45], [100, 100]);
    const chi = S.chiSquareTest([[60, 40], [45, 55]]);
    expect(two.statisticSquared.equals(chi.statistic)).toBe(true);
    expect(two.estimate.equals(Q(3, 20))).toBe(true);
    expect(() => S.proportionTest([0, 0], [5, 5])).toThrow(/pooled variance is 0/);
  });
  test("goodness of fit and independence", () => {
    const g = S.chiSquareTest([22, 17, 21, 18, 12, 30]);
    expect(String(g.statistic)).toBe("91/10");
    expect(g.df).toBe(5);
    expect(mid(g.pValue)).toBeCloseTo(1 - S.chiSquare(5).fast.cdf(9.1), 10);
    const e = S.chiSquareTest([10, 20, 30], { expected: [20, 20, 20] });
    expect(String(e.statistic)).toBe("10");
    const p = S.chiSquareTest([10, 20, 30], { p: ["1/6", "1/3", "1/2"] });
    expect(String(p.statistic)).toBe("0");
    expect(p.pValue.equals(Q(1))).toBe(true);
    const t = S.chiSquareTest([[3, 7, 5], [9, 2, 4]]);
    expect(t.df).toBe(2);
    expect(t.smallExpected).toBe(true);
    // contributions add up to the statistic
    expect(t.contributions.flat().reduce((s, v) => s.add(v), Q(0)).equals(t.statistic)).toBe(true);
    expect(() => S.chiSquareTest([10, 20], { p: ["1/2", "1/3"] })).toThrow(/add up to 5\/6/);
    expect(() => S.chiSquareTest([[1, 0], [2, 0]])).toThrow(/empty/);
  });
  test("correlation: t² from r², the same as the regression slope's t²", () => {
    const xs = [1, 2, 3, 4, 5, 6];
    const ys = [2, 1, 4, 3, 7, 8];
    const r = S.correlationTest(xs, ys);
    const m = S.linearModel(xs, ys);
    expect(r.statisticSquared.equals(m.test("x").statisticSquared)).toBe(true);
    expect(r.rSquared.equals(m.rSquared)).toBe(true);
    const ci = r.interval();
    expect(ci.low.toNumber()).toBeLessThan(mid(r.estimate));
    expect(ci.high.toNumber()).toBeGreaterThan(mid(r.estimate));
    const z = Math.atanh(mid(r.estimate));
    expect(ci.low.toNumber()).toBeCloseTo(Math.tanh(z - 1.959963984540054 / Math.sqrt(3)), 9);
    expect(S.correlationTest([1, 2, 3], [2, 4, 6]).pValue.equals(Q(0))).toBe(true);
  });
});

describe("fitting families", () => {
  const d = [1, 2, 2, 3, 5, 8];
  test("closed forms are exact", () => {
    const n = S.fit(d, "normal");
    expect(String(n.params.mu)).toBe("7/2");
    expect(String(n.params.variance)).toBe("67/12");
    expect(n.distribution).toBe(null);
    expect(S.fit([1, 3], "normal").distribution.toString()).toBe("normal(2, 1)");
    expect(String(S.fit(d, "exponential").params.rate)).toBe("2/7");
    expect(String(S.fit(d, "poisson").params.lambda)).toBe("7/2");
    expect(String(S.fit([0, 1, 1, 0, 1], "bernoulli").params.p)).toBe("3/5");
    expect(String(S.fit([2, 3, 1], "binomial", { trials: 4 }).params.p)).toBe("1/2");
    expect(String(S.fit([1, 2, 3], "geometric").params.p)).toBe("1/2");
    expect(String(S.fit([0, 1, 2], "geometric", { failures: true }).params.p)).toBe("1/2");
    expect(S.fit(d, "uniform").params.b.equals(Q(8))).toBe(true);
    expect(String(S.fit(d, "gamma", { method: "moments" }).params.shape)).toBe("147/67");
    expect(String(S.fit([0.2, 0.5, 0.4, 0.7], "beta", { method: "moments" }).params.a)).toBe("387/130");
  });
  test("the gamma likelihood equation holds at the certified shape", () => {
    const g = S.fit(d, "gamma");
    const k = mid(g.params.shape);
    expect(width(g.params.shape)).toBeLessThan(1e-12);
    const s = Math.log(fmean(d)) - fmean(d.map(Math.log));
    // log k − ψ(k), with ψ from a float series
    const psi = (z) => {
      let r = 0;
      while (z < 8) r -= 1 / z++;
      const z2 = 1 / (z * z);
      return r + Math.log(z) - 0.5 / z - z2 * (1 / 12 - z2 * (1 / 120 - z2 / 252));
    };
    expect(Math.log(k) - psi(k)).toBeCloseTo(s, 8);
    expect(mid(g.params.rate)).toBeCloseTo(k / fmean(d), 12);
  });
  test("lognormal from the logs; refusals", () => {
    const l = S.fit([1, Math.E ** 0 * 2, 4], "lognormal");
    expect(mid(l.params.mu)).toBeCloseTo(Math.log(2), 14);
    expect(() => S.fit(d, "beta")).toThrow(/moments/);
    expect(() => S.fit(d, "weird")).toThrow(/unknown family/);
    expect(() => S.fit([0, 1], "exponential")).toThrow(/positive/);
    expect(() => S.fit([1.5, 2], "poisson")).toThrow(/whole/);
  });
});

describe("regression, exactly", () => {
  const xs = [1, 2, 3, 4, 5, 6, 7, 8];
  const ys = [2.1, 3.9, 6.2, 7.8, 10.1, 12.2, 13.8, 16.1];
  test("the same line as leastSquares, and the normal equations hold exactly", () => {
    const m = S.linearModel(xs, ys);
    const l = S.leastSquares(xs, ys);
    expect(m.coefficient("x").equals(l.slope)).toBe(true);
    expect(m.coefficient("(intercept)").equals(l.intercept)).toBe(true);
    expect(m.rSquared.equals(l.r2)).toBe(true);
    for (const col of m.columns) expect(col.reduce((s, v, i) => s.add(v.multiply(m.residuals[i])), Q(0)).equals(Q(0))).toBe(true);
    expect(m.leverage.reduce((s, v) => s.add(v), Q(0)).equals(Q(2))).toBe(true);
    expect(m.sse.add(m.ssr).equals(m.sst)).toBe(true);
    expect(String(m.predict(10))).toBe("1681/84");
  });
  test("several predictors, recovered exactly from exact data", () => {
    const rows = [[1, 0], [2, 1], [3, 5], [4, 2], [5, 7], [6, 3]];
    const yv = rows.map(([a, b]) => Q(3).add(Q(1, 2).multiply(Q(a))).subtract(Q(2).multiply(Q(b))));
    const m = S.linearModel(rows, yv);
    expect(m.coefficients.map(String)).toEqual(["3", "1/2", "-2"]);
    expect(m.sse.equals(Q(0))).toBe(true);
    expect(m.fTest().pValue.equals(Q(0))).toBe(true);
    const named = S.linearModel({ a: rows.map((r) => r[0]), b: rows.map((r) => r[1]) }, yv);
    expect(named.names).toEqual(["(intercept)", "a", "b"]);
    expect(String(named.predict({ a: 2, b: 2 }))).toBe("0");
    expect(() => S.linearModel([[1, 2], [2, 4], [3, 6], [4, 8]], [1, 2, 3, 5])).toThrow(/x2 is a linear combination of \(intercept\), x1/);
  });
  test("coefficient tests, intervals and prediction intervals hang together", () => {
    const m = S.linearModel(xs, ys);
    const t = m.test("x");
    const ci = m.interval("x");
    expect(m.test("x", { value: ci.low.subtract(Q(1, 1e9)) }).reject()).toBe(true);
    expect(m.test("x", { value: ci.low.add(Q(1, 1e9)) }).reject()).toBe(false);
    expect(t.df).toBe(6);
    const row = m.table()[1];
    expect(row.name).toBe("x");
    expect(mid(row.standardError) ** 2).toBeCloseTo(m.covariance[1][1].toNumber(), 14);
    const mean = m.predictionInterval(10);
    const fresh = m.predictionInterval(10, { kind: "new" });
    expect(fresh.low.compareTo(mean.low)).toBeLessThan(0);
    expect(fresh.high.compareTo(mean.high)).toBeGreaterThan(0);
    // simple regression: F = t² for the slope
    expect(m.fTest().statistic.equals(t.statisticSquared)).toBe(true);
    // Cook's distance from its definition: the change in the fit when a point is left out
    const i = 6;
    const without = S.linearModel(xs.filter((_, j) => j !== i), ys.filter((_, j) => j !== i));
    const change = xs.reduce((s, v, j) => {
      const d = m.fitted[j].subtract(without.predict(v));
      return s.add(d.multiply(d));
    }, Q(0));
    expect(change.divide(Q(2).multiply(m.variance)).equals(m.cooksDistance[i])).toBe(true);
  });
  test("polynomial fits, nested models and dummies", () => {
    const q = S.polyFit([0, 1, 2, 3, 4, 5], [1, 0, 3, 10, 21, 36], 2);
    expect(q.coefficients.map(String)).toEqual(["1", "-3", "2"]);
    expect(String(q.predict(6))).toBe("55");
    const noisy = [1, 0.5, 3, 9.5, 21, 36.5];
    const q2 = S.polyFit([0, 1, 2, 3, 4, 5], noisy, 2);
    const q1 = S.polyFit([0, 1, 2, 3, 4, 5], noisy, 1);
    const c = S.compareModels(q1, q2);
    expect(c.df).toEqual([1, 3]);
    expect(c.statistic.equals(q2.test("x^2").statisticSquared)).toBe(true);
    expect(() => S.compareModels(S.linearModel([0, 1, 2, 3, 4, 5].map((v) => v * v), noisy), q1)).toThrow(/fewer coefficients|not in the span/);
    const g = ["a", "a", "b", "b", "c", "c"];
    expect(S.dummies(g)).toEqual({ b: [0, 0, 1, 1, 0, 0], c: [0, 0, 0, 0, 1, 1] });
    // one-way ANOVA is the F test of the dummies
    const yv = [3, 5, 4, 7, 6, 9];
    const a = S.anova({ a: [3, 5], b: [4, 7], c: [6, 9] });
    expect(S.linearModel(S.dummies(g), yv).fTest().statistic.equals(a.statistic)).toBe(true);
  });
  test("one-way analysis of variance", () => {
    const a = S.anova([[4, 5, 6, 5], [7, 8, 6, 7], [9, 8, 10, 9]]);
    expect(String(a.statistic)).toBe("24");
    expect(a.table.map((r) => r.df)).toEqual([2, 9, 11]);
    expect(String(a.etaSquared)).toBe("16/19");
    expect(mid(a.pValue)).toBeCloseTo(1 - S.fDist(2, 9).fast.cdf(24), 10);
    expect(S.anova([[1, 1], [2, 2]]).pValue.equals(Q(0))).toBe(true);
    expect(() => S.anova([[1, 2]])).toThrow(/two groups/);
  });
});

describe("the procedures do what they claim, by seeded experiment", () => {
  const N = S.normal(0, 1);
  test("the t interval covers the mean 95% of the time", () => {
    const c = S.coverage({ draw: (rng) => sample(S.normal(10, 2), 12, rng), interval: (d) => S.meanInterval(d, { bits: 32 }), truth: 10, runs: 300, seed: 3 });
    expect(c.compare("0.95")).toBe("inside");
  }, 60000);
  test("the t test has its stated size; Welch's test holds it when the variances differ", () => {
    const one = S.rejectionRate({ draw: (rng) => sample(N, 10, rng), test: (d) => S.tTest(d, { bits: 32 }), runs: 300, seed: 5 });
    expect(one.compare("0.05")).toBe("inside");
    const unequal = (rng) => [sample(S.normal(0, 1), 25, rng), sample(S.normal(0, 4), 5, rng)];
    const welch = S.rejectionRate({ draw: unequal, test: ([a, b]) => S.tTest(a, b, { bits: 32 }), runs: 300, seed: 6 });
    expect(welch.compare("0.05")).toBe("inside");
    // the pooled test, given the small sample the large variance, rejects far too often
    const pooled = S.rejectionRate({ draw: unequal, test: ([a, b]) => S.tTest(a, b, { equalVariance: true, bits: 32 }), runs: 300, seed: 6 });
    expect(pooled.interval.low.compareTo(Q(1, 10))).toBeGreaterThan(0);
  }, 120000);
  test("the chi-square test of a fair die has its stated size", () => {
    const counts = (rng) => {
      const c = [0, 0, 0, 0, 0, 0];
      for (let i = 0; i < 60; i++) c[rng.int(6)]++;
      return c;
    };
    const r = S.rejectionRate({ draw: counts, test: (c) => S.chiSquareTest(c, { bits: 32 }), runs: 300, seed: 7 });
    expect(r.compare("0.05")).toBe("inside");
  }, 60000);
  test("a regression slope interval covers the true slope from a seeded model", () => {
    const line = S.model({ truth: [2, 3], x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], noise: S.normal(0, 2), round: "1/100" });
    const c = S.coverage({ model: line, n: 10, interval: (d) => S.linearModel(d.x, d.y).interval("x", { bits: 32 }), truth: 3, runs: 300, seed: 8 });
    expect(c.compare("0.95")).toBe("inside");
  }, 60000);
});
