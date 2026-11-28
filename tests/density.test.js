import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { fromDensity, transform, orderStatistic, uniform, exponential, beta, toNumber, rng } from "../index.js";

const near = (v, x, tol = 1e-6) => expect(Math.abs(toNumber(v) - x)).toBeLessThan(tol);
const q = (s) => new Rational(s);
const eq = (v, s) => expect(String(v)).toBe(String(q(s)));

describe("fromDensity", () => {
  test("normalizes c·x² on [0, 2] (c = 3/8, mean 3/2, variance 3/20)", () => {
    const X = fromDensity("x^2", [0, 2], { normalize: true });
    eq(X.constant, "3/8");
    eq(X.mean, "3/2");
    eq(X.variance, "3/20");
    eq(X.cdf(1), "1/8");
    expect(X.checks).toEqual({ nonnegative: "proven", mass: "exact", massProven: true });
  });
  test("an enclosure of the mass that contains 1 is not reported as a proof", () => {
    // (1 + 10⁻⁵⁰)·cos(x)/sin(1) has mass 1 + 10⁻⁵⁰, inside the certified enclosure of the integral
    const X = fromDensity("(1 + 10^(-50))*cos(x)/sin(1)", [0, 1]);
    expect(X.checks.mass).toBe("enclosure");
    expect(X.checks.massProven).toBe(false);
    expect(X.checks.massEnclosure.low.compareTo(new Rational(1n)) <= 0).toBe(true);
    expect(() => fromDensity("(1 + 10^(-50))*cos(x)/sin(1)", [0, 1], { requireExactMass: true })).toThrow(/does not prove/);
  });
  test("rebuilds uniform(2, 5) exactly", () => {
    const X = fromDensity("1/3", [2, 5]);
    const U = uniform(2, 5);
    eq(X.mean, String(U.mean));
    eq(X.variance, String(U.variance));
    for (const t of ["2", "5/2", "4"]) eq(X.cdf(t), String(U.cdf(t)));
  });
  test("rebuilds beta(2, 3) = 12 x (1 − x)²", () => {
    const X = fromDensity("12*x*(1-x)^2", [0, 1]);
    const B = beta(2, 3);
    eq(X.mean, "2/5");
    eq(X.variance, "1/25");
    for (const t of [0.2, 0.5, 0.9]) near(X.cdf(String(t)), toNumber(B.cdf(String(t))), 1e-9);
  });
  test("rebuilds exponential(1) from exp(-x) with a certified mass", () => {
    const X = fromDensity("exp(-x)", [0, "inf"]);
    const E = exponential(1);
    expect(X.checks.mass).toBe("enclosure");
    near(X.cdf(1), 1 - Math.exp(-1));
    near(X.cdf(3), toNumber(E.cdf(3)));
  });
  test("a piecewise triangle density", () => {
    const T = fromDensity([["x", [0, 1]], ["2 - x", [1, 2]]]);
    eq(T.mean, "1");
    eq(T.variance, "1/6");
    eq(T.cdf("1/2"), "1/8");
    near(T.quantile("1/2"), 1);
  });
  test("refuses a negative formula, a wrong mass and broken pieces", () => {
    expect(() => fromDensity("x - 1", [0, 2])).toThrow(/negative/);
    expect(() => fromDensity("x", [0, 2])).toThrow(/mass is 2/);
    expect(() => fromDensity([["x", [0, 1]], ["1", [2, 3]]])).toThrow(/end to end/);
    expect(() => fromDensity("x", [1, 1])).toThrow(/empty/);
  });
  test("seeded samples match the mean 3/2", () => {
    const X = fromDensity("3/8*x^2", [0, 2]);
    const xs = X.samples(4000, rng(5));
    const m = xs.reduce((a, b) => a + b, 0) / xs.length;
    expect(Math.abs(m - 1.5)).toBeLessThan(0.03);
    expect(xs.every((x) => x >= 0 && x <= 2)).toBe(true);
  });
});

describe("transform", () => {
  test("X² for X uniform on [0, 1]", () => {
    const Y = uniform(0, 1).transform("x^2");
    eq(Y.cdf("1/4"), "1/2");   // g⁻¹(1/4) = 1/2 found exactly
    eq(Y.pdf("1/4"), "1");     // 1/(2·√(1/4))
    near(Y.cdf("1/3"), Math.sqrt(1 / 3), 1e-9);
    eq(Y.mean, "1/3");
    eq(Y.variance, "4/45");
  });
  test("a decreasing g and the function form", () => {
    const Y = transform(uniform(0, 1), "-x");
    eq(Y.cdf("-1/4"), "3/4");
    const Z = transform(exponential(2), "3*x + 1");
    near(Z.cdf(4), 1 - Math.exp(-2));
    near(Z.mean, 2.5);
  });
  test("independent check: P(√X ≤ y) = (y² − 1)/3 for X uniform on [1, 4]", () => {
    const Y = uniform(1, 4).transform("sqrt(x)");
    for (const y of [1.25, 1.5, 1.75]) near(Y.cdf(String(y)), (y * y - 1) / 3, 1e-9);
    expect(() => uniform(0, 4).transform("sqrt(x)")).toThrow(/enclosed/);
  });
  test("refuses non-monotone g and non-affine g on unbounded supports", () => {
    expect(() => uniform(-1, 1).transform("x^2")).toThrow(/monotone/);
    expect(() => exponential(1).transform("x^2")).toThrow(/affine/);
  });
});

describe("orderStatistic", () => {
  test("the max of n uniforms is Beta(n, 1)", () => {
    for (const n of [2, 3, 5]) {
      const M = uniform(0, 1).orderStatistic(n, n);
      eq(M.mean, `${n}/${n + 1}`);
      const B = beta(n, 1);
      for (const t of ["1/4", "3/4"]) near(M.cdf(t), toNumber(B.cdf(t)), 1e-12);
    }
    eq(uniform(0, 1).orderStatistic(3, 3).pdf("1/2"), "3/4");
  });
  test("the min of 3 uniforms is Beta(1, 3), and the median of 3 is Beta(2, 2)", () => {
    eq(orderStatistic(uniform(0, 1), 1, 3).mean, "1/4");
    const Md = orderStatistic(uniform(0, 1), 2, 3);
    eq(Md.mean, "1/2");
    eq(Md.variance, "1/20");
  });
  test("min of n exponentials: P(min ≤ t) = 1 − e^(−nλt)", () => {
    const M = exponential(1).orderStatistic(1, 4);
    near(M.cdf("1/2"), 1 - Math.exp(-2));
  });
  test("seeded simulation of the max of 3 draws from c·x² on [0, 2]", () => {
    const X = fromDensity("3/8*x^2", [0, 2]);
    const M = X.orderStatistic(3, 3);
    const r = rng(9);
    let s = 0;
    for (let k = 0; k < 3000; k++) s += Math.max(...X.samples(3, r));
    expect(Math.abs(s / 3000 - toNumber(M.mean))).toBeLessThan(0.02);
  });
  test("refuses bad k and n", () => {
    expect(() => uniform().orderStatistic(4, 3)).toThrow();
    expect(() => uniform().orderStatistic(0, 3)).toThrow();
  });
});
