import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval } from "@ratmath/core";
import * as S from "../index.js";

const R = (s) => Rational.from(String(s));
const num = S.toNumber;
const width = (v) => (v instanceof RationalInterval ? v.high.subtract(v.low).toNumber() : 0);
const near = (v, x, tol = 1e-11) => Math.abs(num(v) - x) <= tol * Math.max(1, Math.abs(x));
const lo = (v) => (v instanceof RationalInterval ? v.low : v);
const hi = (v) => (v instanceof RationalInterval ? v.high : v);

const ALL = () => [
  S.normal(10, 2), S.lognormal(0, "1/2"), S.gamma(3, 2), S.gamma("1/2", 1), S.chiSquare(5), S.exponential(2), S.weibull("3/2", 2), S.weibull(2, 1),
  S.beta(2, 3), S.beta("1/2", "1/2"), S.studentT(10), S.studentT(1), S.fDist(3, 12), S.uniform(2, 5), S.cauchy(1, 2),
];

describe("values against references", () => {
  test("normal", () => {
    const n = S.normal();
    expect(near(n.cdf("1.96"), 0.9750021048517795)).toBe(true);
    expect(width(n.cdf("1.96"))).toBeLessThan(1e-25);
    expect(near(n.pdf(0), 0.3989422804014327)).toBe(true);
    expect(near(n.quantile("0.975"), 1.959963984540054)).toBe(true);
    expect(String(n.quantile("1/2"))).toBe("0");
    expect(String(n.cdf(0))).toBe("1/2");
    expect(near(n.prob(-1, 1), 0.6826894921370859)).toBe(true);
    expect(near(n.sf(8), 6.220960574271785e-16, 1e-6) || num(n.sf(8)) < 1e-15).toBe(true);
    expect(near(S.normal(10, 2).quantile("0.9"), 12.5631031310892)).toBe(true);
    expect(near(n.quantile("1e-12"), -7.034483825301)).toBe(true);
    expect(String(S.normal(10, 2).mean)).toBe("10");
    expect(String(S.normal(10, 2).variance)).toBe("4");
    expect(String(S.normal(10, 2).sd())).toBe("2");
    expect(String(S.normal(10, 2))).toBe("normal(10, 2)");
    expect(() => S.normal(0, 0)).toThrow(/positive/);
  });
  test("t, chi-square and F", () => {
    expect(near(S.studentT(10).quantile("0.975"), 2.2281388519862744)).toBe(true);
    expect(near(S.studentT(10).pdf(0), 0.3891083839660311)).toBe(true);
    expect(near(S.studentT(1).quantile("0.975"), 12.706204736174694)).toBe(true);
    expect(near(S.studentT(1000).quantile("0.975"), 1.9623390808264078)).toBe(true);
    expect(String(S.studentT(5).cdf(0))).toBe("1/2");
    expect(near(S.studentT(4).cdf(-1), 0.186950483150029)).toBe(true);
    expect(near(S.studentT(4).sf(-1), 1 - 0.186950483150029)).toBe(true);
    expect(near(S.studentT(4).sf(1), 0.186950483150029)).toBe(true);
    expect(String(S.studentT(4).sf(0))).toBe("1/2");
    expect(S.studentT(1).mean).toBe(null);
    expect(S.studentT(2).variance).toBe(null);
    expect(S.studentT(2).sd()).toBe(null);
    expect(String(S.studentT(4).variance)).toBe("2");
    expect(near(S.chiSquare(5).quantile("0.95"), 11.070497693516351)).toBe(true);
    expect(near(S.chiSquare(5).pdf(3), 0.15418032980376925)).toBe(true);
    expect(String(S.chiSquare(2).pdf(0))).toBe("1/2");
    expect(String(S.chiSquare(4).pdf(0))).toBe("0");
    expect(() => S.chiSquare(1).pdf(0)).toThrow(/unbounded/);
    expect(near(S.chiSquare(300).quantile("0.99"), 359.9064259503349)).toBe(true);
    expect(near(S.fDist(3, 12).quantile("0.95"), 3.490294819497605)).toBe(true);
    expect(near(S.fDist(3, 12).pdf(1), S.fDist(3, 12).fast.pdf(1))).toBe(true);
    expect(String(S.fDist(3, 12).mean)).toBe("6/5");
    expect(S.fDist(3, 2).mean).toBe(null);
    expect(S.fDist(3, 4).variance).toBe(null);
    expect(near(S.fDist(3, 12).variance, (2 * 144 * 13) / (3 * 100 * 8))).toBe(true);
    expect(String(S.fDist(2, 5).pdf(0))).toBe("1");
    expect(String(S.fDist(4, 5).pdf(0))).toBe("0");
    expect(() => S.fDist(1, 5).pdf(0)).toThrow(/unbounded/);
  });
  test("gamma, beta, exponential, Weibull, lognormal, uniform, Cauchy", () => {
    expect(near(S.gamma(3, 2).cdf(1), 0.32332358381693654)).toBe(true);
    expect(near(S.gamma(3, 2).quantile("1/2"), 1.3370301568617795)).toBe(true);
    expect(String(S.gamma(1, 3).pdf(0))).toBe("3");
    const b = S.beta(2, 3);
    expect(lo(b.cdf("1/3")).compareTo(R("11/27")) <= 0 && hi(b.cdf("1/3")).compareTo(R("11/27")) >= 0).toBe(true);
    expect(near(b.quantile("1/2"), 0.38572756813238945)).toBe(true);
    expect(near(b.pdf("1/2"), 1.5)).toBe(true);
    expect(String(b.mean)).toBe("2/5");
    expect(String(b.variance)).toBe("1/25");
    expect(String(S.beta(1, 3).pdf(0))).toBe("3");
    expect(String(S.beta(3, 1).pdf(1))).toBe("3");
    expect(String(b.pdf(0))).toBe("0");
    expect(() => S.beta("1/2", 2).pdf(0)).toThrow(/unbounded at 0/);
    expect(() => S.beta(2, "1/2").pdf(1)).toThrow(/unbounded at 1/);
    expect(near(S.exponential(2).quantile("1/2"), Math.LN2 / 2)).toBe(true);
    expect(near(S.exponential(2).cdf(1), 1 - Math.exp(-2))).toBe(true);
    expect(near(S.exponential(2).sf(1), Math.exp(-2))).toBe(true);
    expect(near(S.exponential(2).pdf(1), 2 * Math.exp(-2))).toBe(true);
    expect(String(S.exponential(2).pdf(0))).toBe("2");
    expect(String(S.exponential(2).pdf(-1))).toBe("0");
    const w = S.weibull("3/2", 2);
    expect(near(w.quantile("1/2"), 1.5664395375492327)).toBe(true);
    expect(near(w.cdf(1), 1 - Math.exp(-(0.5 ** 1.5)))).toBe(true);
    expect(near(w.sf(1), Math.exp(-(0.5 ** 1.5)))).toBe(true);
    expect(near(w.mean, 1.8054905859019898)).toBe(true);
    expect(near(w.pdf(1), w.fast.pdf(1))).toBe(true);
    expect(near(S.weibull(2, 1).cdf(1), 1 - Math.exp(-1))).toBe(true);
    expect(String(S.weibull(1, 4).pdf(0))).toBe("1/4");
    expect(String(S.weibull(2, 4).pdf(0))).toBe("0");
    expect(() => S.weibull("1/2", 1).pdf(0)).toThrow(/unbounded/);
    const l = S.lognormal(0, 1);
    expect(near(l.quantile("0.9"), 3.6022244792791563)).toBe(true);
    expect(near(l.cdf(2), 0.7558914042144173)).toBe(true);
    expect(near(l.mean, Math.exp(0.5))).toBe(true);
    expect(near(l.variance, (Math.E - 1) * Math.E)).toBe(true);
    expect(near(l.pdf(1), 0.3989422804014327)).toBe(true);
    expect(near(l.quantile("1/2"), 1)).toBe(true);
    const u = S.uniform(2, 5);
    expect(String(u.quantile("1/3"))).toBe("3");
    expect(String(u.cdf(3))).toBe("1/3");
    expect(String(u.pdf(3))).toBe("1/3");
    expect(String(u.pdf(2))).toBe("1/3");
    expect(String(u.pdf(6))).toBe("0");
    expect(String(u.variance)).toBe("3/4");
    expect(() => S.uniform(1, 1)).toThrow(/a < b/);
    const c = S.cauchy();
    expect(near(c.cdf(1), 0.75)).toBe(true);
    expect(near(c.sf(1), 0.25)).toBe(true);
    expect(near(c.quantile("0.975"), 12.706204736174694)).toBe(true);
    expect(near(c.pdf(0), 1 / Math.PI)).toBe(true);
    expect(c.mean).toBe(null);
  });
});

describe("what is certified", () => {
  test("every quantile enclosure brackets the level", () => {
    for (const d of ALL()) {
      for (const p of ["1/100", "3/10", "1/2", "19/20"]) {
        const q = d.quantile(p, { bits: 40 });
        const level = R(p);
        expect(lo(d.cdf(lo(q), { bits: 80 })).compareTo(level) <= 0).toBe(true);
        expect(hi(d.cdf(hi(q), { bits: 80 })).compareTo(level) >= 0).toBe(true);
        // the solver's promise: width about 2^(−bits) relative to the answer; next to zero (a median of 0) it is absolute
        const size = Math.max(Math.abs(num(q)), 1e-6);
        expect(width(q)).toBeLessThanOrEqual(size * 2 ** -38);
      }
    }
  });
  test("cumulative probabilities enclose, narrow with bits, and agree with the fast side", () => {
    for (const d of ALL()) {
      const x = d.fast.quantile(0.37);
      const q = S.toQ(Number(x.toPrecision(6)));
      const coarse = d.cdf(q, { bits: 30 });
      const fine = d.cdf(q, { bits: 120 });
      expect(lo(coarse).compareTo(lo(fine)) <= 0 && hi(coarse).compareTo(hi(fine)) >= 0).toBe(true);
      expect(width(fine)).toBeLessThan(1e-30);
      expect(Math.abs(num(fine) - d.fast.cdf(q.toNumber()))).toBeLessThan(1e-9);
      expect(Math.abs(num(d.pdf(q)) - d.fast.pdf(q.toNumber()))).toBeLessThan(1e-9 * Math.max(1, d.fast.pdf(q.toNumber())));
      expect(Math.abs(num(d.sf(q)) + num(fine) - 1)).toBeLessThan(1e-20);
      expect(Math.abs(d.fast.quantile(d.fast.cdf(x)) - x)).toBeLessThan(1e-6 * Math.max(1, Math.abs(x)));
    }
  });
  test("the ends of the support, interval arguments and bad arguments", () => {
    const g = S.gamma(3, 2);
    expect(String(g.cdf(-1))).toBe("0");
    expect(String(g.cdf(0))).toBe("0");
    expect(String(g.sf(-1))).toBe("1");
    expect(String(S.uniform(0, 1).cdf(2))).toBe("1");
    expect(String(S.uniform(0, 1).sf(2))).toBe("0");
    expect(String(g.quantile(0))).toBe("0");
    expect(g.quantile(1)).toBe(Infinity);
    expect(S.normal().quantile(0)).toBe(-Infinity);
    expect(String(S.beta(2, 2).quantile(1))).toBe("1");
    const iv = g.cdf("0.9:1.1");
    expect(lo(iv).compareTo(lo(g.cdf("0.9"))) <= 0 && hi(iv).compareTo(hi(g.cdf("1.1"))) >= 0).toBe(true);
    const sv = g.sf("0.9:1.1");
    expect(lo(sv).compareTo(lo(g.sf("1.1"))) <= 0 && hi(sv).compareTo(hi(g.sf("0.9"))) >= 0).toBe(true);
    expect(() => g.prob(2, 1)).toThrow(/a ≤ b/);
    expect(() => g.quantile(2)).toThrow(/between 0 and 1/);
    expect(near(S.exponential(1).quantile("1e-9"), 1.0000000005e-9, 1e-15)).toBe(true);
    expect(near(S.gamma("1/2").quantile("1e-6"), 7.853981634e-13, 1e-18)).toBe(true);
    expect(String(S.normal().median())).toBe("0");
  });
  test("sampling", () => {
    const g = S.rng("continuous");
    for (const d of ALL()) {
      const xs = d.samples(4000, g).sort((a, b) => a - b);
      // the sample median sits near the true one: its rank is binomial(4000, 1/2)
      const below = xs.filter((x) => x < d.fast.quantile(0.5)).length;
      expect(Math.abs(below - 2000)).toBeLessThan(160);
      const tenth = xs.filter((x) => x < d.fast.quantile(0.1)).length;
      expect(Math.abs(tenth - 400)).toBeLessThan(95);
    }
    expect(typeof S.normal().sample(g)).toBe("number");
    expect(S.normal().samples(3, S.rng(5))).toEqual(S.normal().samples(3, S.rng(5)));
    const exact = S.normal().samplesExact(3, S.rng(5));
    expect(exact.every((v) => v instanceof Rational)).toBe(true);
    expect(exact.map(String)).toEqual(S.normal().samplesExact(3, S.rng(5)).map(String));
    expect(S.uniform(0, 1).sampleExact(S.rng(2)) instanceof Rational).toBe(true);
    // an exact draw is the certified quantile of an exact uniform draw
    const u = new Rational((S.rng(8).bits(32) << 1n) | 1n, 1n << 33n);
    const x = S.exponential(1).sampleExact(S.rng(8));
    expect(Math.abs(x.toNumber() + Math.log(1 - u.toNumber()))).toBeLessThan(1e-9);
  });
});
