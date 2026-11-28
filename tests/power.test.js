import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { meanInterval, varianceInterval, power, sampleSize, normal, toNumber, rng, rejectionRate, zTest } from "../index.js";

const near = (v, x, tol = 1e-6) => expect(Math.abs(toNumber(v) - x)).toBeLessThan(tol);

describe("summaries in place of data", () => {
  test("t interval from n, mean, sd (textbook: 23.06 to 27.54)", () => {
    const ci = meanInterval({ n: 16, mean: "25.3", sd: "4.2" });
    near(ci.low, 23.061978, 1e-5);
    near(ci.high, 27.538022, 1e-5);
    expect(ci.method).toBe("t");
  });

  test("a summary gives the same interval as data with that summary (independent check)", () => {
    const xs = [3, 5, 7, 9, 11, 4, 8];
    const n = xs.length, mean = new Rational(47, 7);
    const ss = xs.reduce((s, x) => s.add(new Rational(x).subtract(mean).pow(2)), new Rational(0));
    const variance = ss.divide(new Rational(n - 1));
    const a = meanInterval(xs), b = meanInterval({ n, mean, variance });
    expect(a.low.equals(b.low) && a.high.equals(b.high)).toBe(true);
    const z1 = meanInterval(xs, { sigma: 2 }), z2 = meanInterval({ n, mean }, { sigma: 2 });
    expect(z1.low.equals(z2.low)).toBe(true);
    const v1 = varianceInterval(xs), v2 = varianceInterval({ n, variance });
    expect(v1.low.equals(v2.low) && v1.high.equals(v2.high)).toBe(true);
  });

  test("variance interval from n and s² (textbook: 1.8925 to 13.331)", () => {
    const ci = varianceInterval({ n: 10, variance: 4 });
    near(ci.low, 1.89247, 1e-4);
    near(ci.high, 13.3314, 1e-4);
    expect(varianceInterval({ n: 10, sd: 2 }).low.equals(ci.low)).toBe(true);
  });

  test("refusals", () => {
    expect(() => meanInterval({ n: 16, mean: 3 })).toThrow();
    expect(() => meanInterval({ n: 1, mean: 3, sd: 1 })).toThrow();
    expect(() => varianceInterval({ n: 10, sd: 2, variance: 4 })).toThrow();
  });
});

describe("power", () => {
  test("textbook z test: 0.6388", () => {
    near(power({ mu0: 100, mu: 106, sigma: 15, n: 25, alternative: "greater" }).power, 0.63876, 1e-5);
  });

  test("matches 1 − Φ(z₁₋α − δ) computed in floats (independent check)", () => {
    const Phi = (x) => toNumber(normal().cdf(String(x)));
    for (const [mu, n] of [[103, 16], [106, 25], [110, 9]]) {
      const delta = ((mu - 100) * Math.sqrt(n)) / 15;
      near(power({ mu0: 100, mu, sigma: 15, n, alternative: "greater" }).power, 1 - Phi((1.6448536269514722 - delta).toFixed(15)), 1e-9);
      near(power({ mu0: 100, mu: 200 - mu, sigma: 15, n, alternative: "less" }).power, 1 - Phi((1.6448536269514722 - delta).toFixed(15)), 1e-9);
    }
  });

  test("two-sided power at the null is alpha", () => {
    near(power({ mu0: 0, mu: 0, sigma: 1, n: 10 }).power, 0.05, 1e-12);
  });

  test("with no effect, n = 1 reaches a target power up to alpha and nothing reaches more", () => {
    expect(sampleSize({ power: "0.04", mu0: 0, mu: 0, sigma: 1 }).n).toBe(1);
    expect(sampleSize({ power: "0.04", test: "proportion", p0: "0.3", p: "0.3" }).n).toBe(1);
    expect(() => sampleSize({ power: "0.5", mu0: 0, mu: 0, sigma: 1 })).toThrow(/never reaches/);
  });

  test("a sample size above 2^53 comes back as an exact bigint", () => {
    // (1.959963…/10⁻⁹)²: both ends of the enclosure have the same ceiling
    const r = sampleSize({ margin: "1/1000000000", sigma: 1 });
    expect(r.n).toBe(3841458820694125959n);
    expect(sampleSize({ margin: 2, sigma: 10 }).n).toBe(97);
  });

  test("agrees with simulation (seeded)", () => {
    const r = rejectionRate({
      draw: (rand) => normal(106, 15).samples(25, rand).map((x) => x.toFixed(3)),
      test: (xs) => zTest(xs, { mu: 100, sigma: 15, alternative: "greater" }),
      alpha: "0.05", runs: 400, seed: 11,
    });
    const p = toNumber(power({ mu0: 100, mu: 106, sigma: 15, n: 25, alternative: "greater" }).power);
    // the exact interval for the simulated rate contains the closed-form power
    expect(r.compare(String(p.toFixed(6)))).toBe("inside");
  });

  test("proportion test by the normal approximation", () => {
    near(power({ test: "proportion", p0: "0.5", p: "0.6", n: 100, alternative: "greater" }).power, 0.6415, 1e-3);
    expect(() => power({ test: "t", n: 5 })).toThrow();
  });
});

describe("sampleSize", () => {
  test("margins of error", () => {
    expect(sampleSize({ margin: "0.03" }).n).toBe(1068);
    expect(sampleSize({ margin: "0.05", p: "0.2" }).n).toBe(246); // 1.96²·0.16/0.0025 = 245.86
    expect(sampleSize({ margin: 2, sigma: 10 }).n).toBe(97);
    expect(sampleSize({ margin: "0.03", level: "0.99" }).n).toBe(1844); // 2.5758²/0.0036 = 1843.0
  });

  test("power targets: the smallest n, checked on both sides", () => {
    const spec = { mu0: 100, mu: 106, sigma: 15, alternative: "greater" };
    const s = sampleSize({ power: "0.8", ...spec });
    expect(s.n).toBe(39); // ((1.6449 + 0.8416)·15/6)² = 38.64
    expect(toNumber(power({ ...spec, n: 39 }).power)).toBeGreaterThanOrEqual(0.8);
    expect(toNumber(power({ ...spec, n: 38 }).power)).toBeLessThan(0.8);
    expect(sampleSize({ power: "0.9", mu0: 0, mu: "0.5", sigma: 1 }).n).toBe(43);
    const sp = sampleSize({ power: "0.8", test: "proportion", p0: "0.5", p: "0.6", alternative: "greater" });
    expect(toNumber(power({ test: "proportion", p0: "0.5", p: "0.6", alternative: "greater", n: sp.n }).power)).toBeGreaterThanOrEqual(0.8);
    expect(toNumber(power({ test: "proportion", p0: "0.5", p: "0.6", alternative: "greater", n: sp.n - 1 }).power)).toBeLessThan(0.8);
    expect(() => sampleSize({ power: "0.8", mu0: 1, mu: 1, sigma: 1 })).toThrow();
  });
});
