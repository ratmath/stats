import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval, Integer } from "@ratmath/core";
import { permutations, combinations } from "@ratmath/discrete";
import * as S from "../index.js";

const R = (s) => Rational.from(String(s));
const str = (v) => String(v);

describe("reading numbers", () => {
  test("JS numbers are read by their decimal form", () => {
    expect(str(S.toQ(0.1))).toBe("1/10");
    expect(str(S.toQ(-2.5))).toBe("-5/2");
    expect(str(S.toQ(1e-7))).toBe("1/10000000");
    expect(str(S.toQ(1.5e21))).toBe("1500000000000000000000");
    expect(str(S.toQ("1e-3"))).toBe("1/1000");
    expect(str(S.toQ("2.5E2"))).toBe("250");
    expect(str(S.toQ(7n))).toBe("7");
    expect(str(S.toQ("3/4"))).toBe("3/4");
    expect(() => S.toQ(NaN)).toThrow(/finite/);
    expect(() => S.toQ("abc")).toThrow(/cannot read/);
    expect(() => S.toQ({})).toThrow(/must be a number/);
    expect(() => S.toQ(null)).toThrow(/must be a number/);
  });
  test("data values may be intervals", () => {
    expect(S.toValue("1.2:1.4") instanceof RationalInterval).toBe(true);
    expect(str(S.toValue([1, 2]))).toBe("1:2");
    expect(str(S.toValue([2, 1]))).toBe("1:2");
    expect(str(S.toValue(new RationalInterval(R(3), R(3))))).toBe("3");
    expect(str(S.toValue("1.5[+-0.1]"))).toBe("149/100:151/100");
    expect(() => S.toValue("1:x")).toThrow(/cannot read/);
  });
  test("show, compare and sqrt", () => {
    expect(S.show(R("1/3"), 4)).toBe("0.3333");
    expect(S.show(R("-2/3"), 3)).toBe("-0.667");
    expect(S.show(R("5"))).toBe("5");
    expect(S.show(new RationalInterval(R("0.12345"), R("0.12349")), 6)).toBe("0.1235…");
    expect(S.show(new RationalInterval(R("0.1"), R("0.9")), 2)).toBe("[0.10, 0.90]");
    expect(S.show(new RationalInterval(R("0.123456701"), R("0.123456702")), 6)).toBe("0.123457");
    expect(S.show(0.1234567891, 4)).toBe("0.12346");
    expect(S.show("text")).toBe("text");
    expect(S.compare(R("0.03"), "0.05")).toBe("below");
    expect(S.compare(R("0.05"), "0.05")).toBe("equal");
    expect(S.compare(new RationalInterval(R("0.04"), R("0.06")), "0.05")).toBe("too close to call");
    expect(S.compare(new RationalInterval(R("0.06"), R("0.07")), "0.05")).toBe("above");
    expect(str(S.sqrt(R("9/4")))).toBe("3/2");
    const r = S.sqrt(R(2));
    expect(r.low.multiply(r.low).compareTo(R(2)) < 0 && r.high.multiply(r.high).compareTo(R(2)) > 0).toBe(true);
    expect(() => S.sqrt(R(-1))).toThrow(/negative/);
    expect(S.toNumber(new RationalInterval(R(1), R(2)))).toBe(1.5);
  });
});

describe("finite probability", () => {
  const two = S.dice(2);
  test("events, conditioning and independence", () => {
    expect(two.size).toBe(36);
    expect(str(two.prob(([a, b]) => a + b === 7))).toBe("1/6");
    expect(str(two.prob(([a, b]) => a + b >= 10, ([a]) => a === 6))).toBe("1/2");
    expect(str(two.p([3, 4]))).toBe("1/36");
    expect(str(two.p([7, 7]))).toBe("0");
    expect(str(two.prob(new Set([[1, 1], [6, 6]])))).toBe("1/18");
    expect(str(two.prob([1, 1]))).toBe("1/36");
    expect(two.independent(([a]) => a % 2 === 0, ([a, b]) => a + b === 7)).toBe(true);
    expect(two.independent(([a]) => a, ([, b]) => b)).toBe(true);
    expect(two.independent(([a]) => a, ([a, b]) => a + b)).toBe(false);
    expect(two.given(([a, b]) => a === b).size).toBe(6);
    expect(two.filter(([a]) => a === 1).size).toBe(6);
    expect(() => two.given(() => false)).toThrow(/probability 0/);
  });
  test("random variables, moments and quantiles", () => {
    const sum = two.map(([a, b]) => a + b);
    expect(str(sum.mean())).toBe("7");
    expect(str(sum.expect((s) => s * s))).toBe("329/6");
    expect(str(sum.variance())).toBe("35/6");
    expect(str(sum.moment(3, { central: true }))).toBe("0");
    expect(str(sum.cdf(4))).toBe("1/6");
    expect(str(sum.quantile("1/2"))).toBe("7");
    expect(str(sum.median())).toBe("7");
    expect(sum.modes()).toEqual([7]);
    expect(sum.support.map(([x]) => str(x)).join(" ")).toBe("2 3 4 5 6 7 8 9 10 11 12");
    const sd = sum.sd();
    expect(sd instanceof RationalInterval && Math.abs(S.toNumber(sd) - Math.sqrt(35 / 6)) < 1e-12).toBe(true);
    expect(str(S.finite([1, 3]).sd())).toBe("1");
    expect(str(two.covariance(([a]) => a, ([a, b]) => a + b))).toBe("35/12");
    expect(Math.abs(S.toNumber(two.correlation(([a]) => a, ([a, b]) => a + b)) - Math.SQRT1_2)).toBeLessThan(1e-12);
    expect(str(two.correlation(([a]) => a, ([a]) => -a))).toBe("-1");
    expect(() => two.correlation(([a]) => a, () => 1)).toThrow(/constant/);
    expect(() => S.coin().mean()).toThrow(/needs numbers/);
    expect(str(two.marginal(0).p(3))).toBe("1/6");
    expect(str(two.joint(([a]) => a, ([a, b]) => a * b).p([2, 6]))).toBe("1/36");
    expect(() => sum.marginal(0)).toThrow(/tuples/);
  });
  test("sums, products and powers", () => {
    expect(str(S.die().sumOf(10).p(35))).toBe("7631/104976");
    expect(str(S.die().plus(S.die()).p(7))).toBe("1/6");
    // Core's Integer (what R`1` gives) counts as a number
    expect(str(S.die().mean((d) => (d % 2 === 0 ? new Integer(1n) : new Rational(1n, 2n))))).toBe("3/4");
    expect(str(S.die().map((d) => new Integer(BigInt(d % 2))).p(new Integer(1n)))).toBe("1/2");
    const heads = S.coin("1/3").power(3).map((t) => t.filter((c) => c === "H").length);
    expect(str(heads.p(2))).toBe("2/9");
    expect(S.coin().product(S.die(), ["x", "y"]).size).toBe(24);
    expect(S.dice(9).size).toBe(10077696);
    expect(() => S.dice(9).entries).toThrow(/limit/);
    expect(str(heads)).toContain("3: 1/27");
    expect(S.dice(2).toString()).toContain("36 outcomes");
  });
  test("building distributions", () => {
    expect(str(S.finite({ H: "1/2", T: "1/2" }).p("H"))).toBe("1/2");
    expect(str(S.finite(["a", "a", "b"]).p("a"))).toBe("2/3");
    expect(str(S.finite(["a", "b"], [1, 3], { normalize: true }).p("b"))).toBe("3/4");
    expect(str(S.finite(new Map([[1, "1/4"], [2, "3/4"]])).mean())).toBe("7/4");
    expect(() => S.finite(["a", "b"], ["1/2", "1/3"])).toThrow(/add up to 5\/6/);
    expect(() => S.finite(["a"], [1, 2])).toThrow(/one weight/);
    expect(() => S.finite([])).toThrow(/at least one/);
    expect(() => S.finite(5)).toThrow(/needs an array/);
    expect(() => S.finite({ a: -1, b: 2 })).toThrow(/negative/);
    expect(() => S.finite({ a: 0 })).toThrow(/positive probability/);
    expect(str(S.finite({ a: 2, b: 2 }, { normalize: true }).p("a"))).toBe("1/2");
    // 0.5 and 1/2 are the same outcome
    expect(S.finite([0.5, R("1/2"), 1]).size).toBe(2);
    const f = S.finite(S.finite([1, 2]));
    expect(f.size).toBe(2);
    // derangements of 4 by enumeration
    expect(str(S.finite(permutations(4)).prob((p) => p.every((v, i) => v !== i + 1)))).toBe("3/8");
    expect(() => S.finite(combinations(40, 20))).toThrow(/limit/);
    expect(S.deck().length).toBe(52);
    expect(str(S.draw(S.deck(), 2).prob((h) => h.every((c) => c.startsWith("A"))))).toBe("1/221");
    expect(S.draw(["a", "b", "c"], 2, { ordered: true }).size).toBe(6);
    expect(S.draw(["a", "b", "c"], 2, { ordered: true, replace: true }).size).toBe(9);
    expect(() => S.draw(["a", "b"], 2, { replace: true })).toThrow(/not equally likely/);
    expect(S.finite([1, 2]).table()).toHaveLength(2);
    expect(S.finite([3, 1, 2]).table({ sort: true }).map((r) => r.outcome)).toEqual([1, 2, 3]);
    expect(S.finite(["b", "a"]).table({ sort: (p, q) => p.localeCompare(q) })[0].outcome).toBe("a");
  });
  test("Bayes' rule", () => {
    const b = S.bayes({ sick: "1/100", well: "99/100" }, { sick: "95/100", well: "5/100" });
    expect(str(b.posterior.p("sick"))).toBe("19/118");
    expect(str(b.evidence)).toBe("59/1000");
    expect(b.table).toHaveLength(2);
    const c = S.bayes(S.finite([1, 2, 3]), (h) => R(h).divide(R(3)));
    expect(str(c.posterior.p(3))).toBe("1/2");
    const d = S.bayes(S.finite(["a", "b"]), new Map([["a", 1], ["b", 0]]));
    expect(str(d.posterior.p("a"))).toBe("1");
    expect(() => S.bayes({ a: 1 }, { a: 0 })).toThrow(/probability 0/);
    expect(() => S.bayes({ a: 1 }, {})).toThrow(/no likelihood/);
  });
  test("sampling follows the exact probabilities", () => {
    const g = S.rng(7);
    const d = S.finite({ a: "1/10", b: "9/10" });
    const draws = d.samples(4000, g);
    const a = draws.filter((v) => v === "a").length;
    expect(a).toBeGreaterThan(330);
    expect(a).toBeLessThan(470);
    expect(d.samples(5, S.rng(1))).toEqual(d.samples(5, S.rng(1)));
  });
});

describe("outcome identity", () => {
  test("tuples whose members contain the separator stay apart (Astra review)", () => {
    const d = S.finite([["a,sb", "c"], ["a", "b,sc"]]);
    expect(d.entries.length).toBe(2);
    expect(str(d.p(["a,sb", "c"]))).toBe("1/2");
    expect(str(d.p(["a", "b,sc"]))).toBe("1/2");
    expect(str(d.prob((o) => o[0] === "a"))).toBe("1/2");
  });
  test("sets, nested tuples and objects keep distinct members apart", () => {
    expect(S.finite([new Set(["a,sb"]), new Set(["a", "b"])]).entries.length).toBe(2);
    expect(S.finite([[["a"], "b"], [["a", "b"]]]).entries.length).toBe(2);
    expect(S.finite([{ x: 1n }, { x: "1n" }]).entries.length).toBe(2);
    expect(S.finite([{ x: R("1/2") }, { x: "1/2" }]).entries.length).toBe(2);
  });
  test("equal numbers still agree whatever their type", () => {
    expect(S.finite([[1, 2], [1n, new Integer(2n)], [R(1), 2]]).entries.length).toBe(1);
    expect(S.finite([new Set([1, 2]), new Set([2n, 1])]).entries.length).toBe(1);
    expect(S.finite([{ x: 5 }, { x: 5n }]).entries.length).toBe(1);
  });
});

describe("discrete distributions", () => {
  test("binomial", () => {
    const b = S.binomial(10, "3/10");
    expect(str(b.pmf(3))).toBe("66706983/250000000");
    expect(str(b.pmf(11))).toBe("0");
    expect(str(b.pmf(2.5))).toBe("0");
    expect(str(b.cdf(10))).toBe("1");
    expect(str(b.cdf(-1))).toBe("0");
    expect(str(b.cdf(2).add(b.sf(2)))).toBe("1");
    expect(str(b.prob(2, 4))).toBe(str(b.pmf(2).add(b.pmf(3)).add(b.pmf(4))));
    expect(str(b.mean)).toBe("3");
    expect(str(b.variance)).toBe("21/10");
    expect(b.mode).toBe(3);
    expect(b.quantile("1/2")).toBe(3);
    expect(b.median()).toBe(3);
    expect(b.quantile(1)).toBe(10);
    expect(str(b.toFinite().mean())).toBe("3");
    expect(String(b)).toBe("binomial(10, 3/10)");
    expect(Math.abs(S.toNumber(b.sd()) - Math.sqrt(2.1))).toBeLessThan(1e-12);
    expect(() => S.binomial(-1, "1/2")).toThrow(/whole number/);
    expect(() => S.binomial(3, 2)).toThrow(/between 0 and 1/);
  });
  test("the others, against their definitions", () => {
    expect(str(S.bernoulli("1/4").pmf(1))).toBe("1/4");
    expect(S.bernoulli("3/4").mode).toBe(1);
    const g = S.geometric("1/6");
    expect(str(g.pmf(3))).toBe("25/216");
    expect(str(g.cdf(3))).toBe("91/216");
    expect(str(g.mean)).toBe("6");
    expect(g.quantile("1/2")).toBe(4);
    expect(g.quantile(1)).toBe(Infinity);
    expect(() => g.toFinite()).toThrow(/infinitely many/);
    const gf = S.geometric("1/6", { failures: true });
    expect(str(gf.pmf(2))).toBe("25/216");
    expect(str(gf.mean)).toBe("5");
    expect(() => S.geometric(0)).toThrow(/p > 0/);
    const nb = S.negativeBinomial(3, "1/2");
    expect(str(nb.pmf(2))).toBe("3/16");
    expect(str(nb.mean)).toBe("3");
    expect(str(nb.cdf(1))).toBe("5/16");
    expect(nb.mode).toBe(2);
    expect(S.negativeBinomial(1, "1/2").mode).toBe(0);
    const h = S.hypergeometric(52, 4, 5);
    expect(str(h.pmf(2))).toBe("2162/54145");
    expect(str(h.mean)).toBe("5/13");
    expect(h.min).toBe(0);
    expect(h.max).toBe(4);
    expect(str(h.toFinite().variance())).toBe(str(h.variance));
    expect(() => S.hypergeometric(5, 6, 2)).toThrow(/cannot exceed/);
    const u = S.discreteUniform(-2, 3);
    expect(str(u.pmf(0))).toBe("1/6");
    expect(str(u.cdf(0))).toBe("1/2");
    expect(str(u.mean)).toBe("1/2");
    expect(str(u.variance)).toBe("35/12");
    expect(() => S.discreteUniform(3, 1)).toThrow(/a ≤ b/);
  });
  test("Poisson probabilities are certified", () => {
    const p = S.poisson(4);
    expect(p.kind).toBe("certified");
    const v = p.pmf(3);
    expect(v instanceof RationalInterval).toBe(true);
    expect(Math.abs(S.toNumber(v) - 0.19536681481316454)).toBeLessThan(1e-15);
    expect(v.high.subtract(v.low).toNumber()).toBeLessThan(1e-25);
    expect(Math.abs(S.toNumber(p.cdf(6)) - 0.8893260215974266)).toBeLessThan(1e-15);
    expect(S.toNumber(p.pmf(3, { bits: 20 })) > 0.19).toBe(true);
    expect(p.quantile("0.95")).toBe(8);
    expect(str(p.mean)).toBe("4");
    expect(() => p.toFinite()).toThrow(/infinitely many/);
    expect(() => S.poisson(0)).toThrow(/λ > 0/);
    expect(S.toNumber(p.prob(0, 2))).toBeCloseTo(0.2381033, 6);
  });
  test("samplers agree with the probabilities", () => {
    const g = S.rng("samplers");
    const check = (d, n = 4000) => {
      const xs = d.samples(n, g);
      const m = xs.reduce((a, b) => a + b, 0) / n;
      const se = Math.sqrt(d.variance.toNumber() / n);
      expect(Math.abs(m - d.mean.toNumber())).toBeLessThan(5 * se);
    };
    check(S.binomial(12, "1/3"));
    check(S.bernoulli("1/5"));
    check(S.geometric("1/4"));
    check(S.geometric("1/4", { failures: true }));
    check(S.negativeBinomial(3, "2/5"));
    check(S.hypergeometric(30, 12, 7));
    check(S.discreteUniform(1, 6));
    check(S.poisson("7/2"), 2000);
    expect(S.poisson(3).samples(5, S.rng(9))).toEqual(S.poisson(3).samples(5, S.rng(9)));
  });
});

describe("descriptive statistics", () => {
  const data = [2, 4, 4, 4, 5, 5, 7, 9];
  test("centre and spread, exactly", () => {
    expect(S.count(data)).toBe(8);
    expect(str(S.sum(data))).toBe("40");
    expect(str(S.mean(data))).toBe("5");
    expect(str(S.variance(data, { population: true }))).toBe("4");
    expect(str(S.sd(data, { population: true }))).toBe("2");
    expect(str(S.variance(data))).toBe("32/7");
    expect(S.sd(data) instanceof RationalInterval).toBe(true);
    expect(str(S.mean([0.1, 0.2, 0.3]))).toBe("1/5");
    expect(str(S.weightedMean([1, 2, 3], [1, 1, 2]))).toBe("9/4");
    expect(str(S.moment(data, 2))).toBe("29");
    expect(str(S.kurtosis(data))).toBe("89/32");
    expect(str(S.kurtosis(data, { excess: true }))).toBe("-7/32");
    expect(Math.abs(S.toNumber(S.skewness(data)) - 0.65625)).toBeLessThan(1e-12);
    expect(S.toNumber(S.skewness([1, 2, 3, 4, 10])) > 0).toBe(true);
    expect(S.toNumber(S.skewness([-10, 1, 2, 3, 4])) < 0).toBe(true);
    expect(str(S.min(data))).toBe("2");
    expect(str(S.max(data))).toBe("9");
    expect(str(S.range(data))).toBe("7");
    expect(str(S.mad(data))).toBe("1/2");
    expect(str(S.trimmedMean(data, "1/8"))).toBe("29/6");
    expect(S.modes(data).map(str)).toEqual(["4"]);
    expect(S.modes([1, 1, 2, 2, 3]).map(str)).toEqual(["1", "2"]);
    expect(() => S.mean([])).toThrow(/empty/);
    expect(() => S.mean("abc")).toThrow(/list of numbers/);
    expect(() => S.variance([1])).toThrow(/at least two/);
    expect(() => S.skewness([1, 1])).toThrow(/equal/);
    expect(() => S.kurtosis([1, 1])).toThrow(/equal/);
    expect(() => S.trimmedMean(data, "1/2")).toThrow(/below 1\/2/);
    expect(() => S.weightedMean([1, 2], [1])).toThrow(/one weight/);
  });
  test("textbook quartiles: medians of the halves", () => {
    const odd = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    expect(S.quartiles(odd, { method: "halves" }).map(str)).toEqual(["5/2", "5", "15/2"]);
    expect(S.quartiles(odd, { method: "tukey" }).map(str)).toEqual(["3", "5", "7"]);
    const even = [3, 7, 8, 12, 13, 14, 18, 21];
    expect(S.quartiles(even, { method: "halves" }).map(str)).toEqual(["15/2", "25/2", "16"]);
    expect(str(S.iqr(even, { method: "halves" }))).toBe("17/2");
    expect(S.fiveNumber([6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49], { method: "halves" }).map(str)).toEqual(["6", "15", "40", "43", "49"]);
    expect(S.fiveNumber([6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49], { method: "tukey" }).map(str)).toEqual(["6", "51/2", "40", "85/2", "49"]);
    expect(str(S.quartiles([4], { method: "halves" })[0])).toBe("4");
    expect(() => S.quantile(odd, "1/10", { method: "halves" })).toThrow(/only the quartiles/);
  });
  test("quantiles by every method", () => {
    const xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    expect(str(S.median(xs))).toBe("11/2");
    expect(str(S.quantile(xs, "1/4"))).toBe("13/4");
    expect(str(S.quantile(xs, "1/4", { method: "weibull" }))).toBe("11/4");
    expect(str(S.quantile(xs, "1/4", { method: "hazen" }))).toBe("3");
    expect(str(S.quantile(xs, "1/4", { method: "median-unbiased" }))).toBe("35/12");
    expect(str(S.quantile(xs, "1/4", { method: "inverse" }))).toBe("3");
    expect(str(S.quantile(xs, "3/10", { method: "inverse" }))).toBe("3");
    expect(str(S.quantile(xs, "1/4", { method: "lower" }))).toBe("3");
    expect(str(S.quantile(xs, "1/4", { method: "higher" }))).toBe("4");
    expect(str(S.quantile(xs, "1/4", { method: "nearest" }))).toBe("3");
    expect(str(S.quantile(xs, "1/3", { method: "nearest" }))).toBe("4");
    expect(str(S.quantile(xs, "1/4", { method: "midpoint" }))).toBe("7/2");
    expect(str(S.quantile(xs, "1/3", { method: "midpoint" }))).toBe("4");
    expect(str(S.quantile(xs, 0))).toBe("1");
    expect(str(S.quantile(xs, 1))).toBe("10");
    expect(str(S.quantile(xs, 0, { method: "weibull" }))).toBe("1");
    expect(str(S.quantile(xs, 1, { method: "weibull" }))).toBe("10");
    expect(str(S.quantile(xs, 0, { method: "inverse" }))).toBe("1");
    expect(S.quartiles(xs).map(str)).toEqual(["13/4", "11/2", "31/4"]);
    expect(str(S.iqr(xs))).toBe("9/2");
    expect(S.fiveNumber(xs).map(str)).toEqual(["1", "13/4", "11/2", "31/4", "10"]);
    expect(() => S.quantile(xs, "1/2", { method: "magic" })).toThrow(/unknown quantile method/);
    expect(() => S.quantile(xs, 2)).toThrow(/between 0 and 1/);
  });
  test("tables and histograms", () => {
    const f = S.frequencies(data);
    expect(f.map((r) => r.count)).toEqual([1, 3, 2, 1, 1]);
    expect(str(f[1].proportion)).toBe("3/8");
    expect(str(f[4].cumulative)).toBe("1");
    const h = S.histogram(data, { width: 2, start: 2 });
    expect(h.map((r) => r.count)).toEqual([1, 5, 1, 1]);
    expect(str(h[1].density)).toBe("5/16");
    expect(S.histogram(data, { edges: [0, 5, 10] }).map((r) => r.count)).toEqual([4, 4]);
    expect(S.histogram(data, { edges: [0, 5] }).outside).toBe(2);
    expect(S.histogram(data, { bins: 7 }).length).toBe(7);
    expect(S.histogram(data).length).toBe(4);
    expect(S.histogram([3, 3, 3]).length).toBe(1);
    expect(() => S.histogram(data, { edges: [1] })).toThrow(/two edges/);
    expect(() => S.histogram(data, { edges: [1, 1] })).toThrow(/increase/);
    expect(() => S.histogram(data, { width: 0 })).toThrow(/positive/);
  });
  test("two variables", () => {
    const x = [1, 2, 3, 4, 5];
    const y = [2, 4, 5, 4, 5];
    expect(str(S.covariance(x, y))).toBe("3/2");
    expect(str(S.covariance(x, y, { population: true }))).toBe("6/5");
    const r = S.correlation(x, y);
    expect(Math.abs(S.toNumber(r) - 0.7745966692414834)).toBeLessThan(1e-13);
    expect(str(S.correlation(x, [2, 4, 6, 8, 10]))).toBe("1");
    expect(str(S.correlation(x, [5, 4, 3, 2, 1]))).toBe("-1");
    expect(S.toNumber(S.correlation(x, [5, 4, 5, 2, 1])) < 0).toBe(true);
    const fit = S.leastSquares(x, y);
    expect(str(fit.slope)).toBe("3/5");
    expect(str(fit.intercept)).toBe("11/5");
    expect(str(fit.r2)).toBe("3/5");
    expect(str(fit.predict(6))).toBe("29/5");
    expect(str(fit.sse)).toBe("12/5");
    expect(str(fit.residuals.reduce((s, v) => s.add(v), R(0)))).toBe("0");
    expect(str(S.leastSquares([1, 2], [3, 3]).r2)).toBe("1");
    expect(() => S.leastSquares([1, 1], [1, 2])).toThrow(/no line/);
    expect(() => S.leastSquares([1], [1])).toThrow(/two points/);
    expect(() => S.covariance([1, 2], [1])).toThrow(/x has 2/);
    expect(() => S.covariance([1], [1])).toThrow(/two pairs/);
    expect(() => S.correlation([1, 1], [1, 2])).toThrow(/constant/);
  });
  test("summary", () => {
    const s = S.summary([1, 2, 3, 4, 100]);
    expect(s.n).toBe(5);
    expect(str(s.median)).toBe("3");
    expect(s.outliers.map(str)).toEqual(["100"]);
    expect(str(s.iqr)).toBe("2");
    expect(S.summary([5]).variance).toBeUndefined();
  });
  test("interval data give enclosures", () => {
    const xs = ["1.0:1.2", "1.9:2.1", [2.8, 3.2]];
    const m = S.mean(xs);
    expect(str(m)).toBe("19/10:13/6");
    expect(str(S.median(xs))).toBe("19/10:21/10");
    expect(str(S.min(xs))).toBe("1:6/5");
    expect(str(S.max(xs))).toBe("14/5:16/5");
    expect(str(S.range(xs))).toBe("8/5:11/5");
    const v = S.variance(xs);
    // every choice of points gives a variance inside the enclosure
    for (const pick of [[1, 1.9, 2.8], [1.2, 2.1, 3.2], [1.1, 2, 3], [1.2, 1.9, 3.2]]) {
      const exact = S.variance(pick);
      expect(v.low.compareTo(exact) <= 0 && v.high.compareTo(exact) >= 0).toBe(true);
    }
    expect(S.sd(xs) instanceof RationalInterval).toBe(true);
    expect(S.summary(xs).outliers).toBeUndefined();
    expect(str(S.trimmedMean(xs, 0))).toBe(str(m));
    expect(() => S.frequencies(xs)).toThrow(/interval data/);
    expect(() => S.moment(xs, 2)).toThrow(/interval data/);
    expect(str(S.range(["1:3", "2:4"]))).toBe("0:3");
    expect(str(S.iqr(xs)).includes(":")).toBe(true);
  });
  test("the fast versions agree", () => {
    const x = [1.5, 2.25, 3.1, 4.8, 5.05, 7.7];
    const y = [2.1, 2.9, 4.2, 5.1, 5.8, 9.4];
    expect(S.fast.mean(x)).toBeCloseTo(S.mean(x).toNumber(), 12);
    expect(S.fast.sum(x)).toBeCloseTo(S.sum(x).toNumber(), 12);
    expect(S.fast.variance(x)).toBeCloseTo(S.variance(x).toNumber(), 12);
    expect(S.fast.variance(x, { population: true })).toBeCloseTo(S.variance(x, { population: true }).toNumber(), 12);
    expect(S.fast.sd(x)).toBeCloseTo(S.toNumber(S.sd(x)), 12);
    expect(S.fast.median(x)).toBeCloseTo(S.median(x).toNumber(), 12);
    expect(S.fast.quantile(x, 0.3)).toBeCloseTo(S.quantile(x, "0.3").toNumber(), 12);
    expect(S.fast.quantile(x, 1)).toBe(7.7);
    expect(S.fast.min(x)).toBe(1.5);
    expect(S.fast.max(x)).toBe(7.7);
    expect(S.fast.trimmedMean(x, 0.2)).toBeCloseTo(S.trimmedMean(x, "0.2").toNumber(), 12);
    expect(S.fast.covariance(x, y)).toBeCloseTo(S.covariance(x, y).toNumber(), 12);
    expect(S.fast.correlation(x, y)).toBeCloseTo(S.toNumber(S.correlation(x, y)), 12);
    const f = S.fast.leastSquares(x, y);
    const e = S.leastSquares(x, y);
    expect(f.slope).toBeCloseTo(e.slope.toNumber(), 12);
    expect(f.intercept).toBeCloseTo(e.intercept.toNumber(), 12);
    expect(f.r2).toBeCloseTo(e.r2.toNumber(), 12);
    expect(f.predict(2)).toBeCloseTo(e.predict(2).toNumber(), 12);
    expect(S.fast.leastSquares([1, 2], [3, 3]).r2).toBe(1);
  });
});

describe("exact tests", () => {
  test("the binomial test", () => {
    const t = S.binomialTest(7, 20);
    expect(str(t.pValue)).toBe("34495/131072");
    expect(t.kind).toBe("exact");
    expect(t.reject()).toBe(false);
    expect(t.decision("0.3")).toBe("reject");
    expect(String(t)).toContain("exact binomial test: p = 0.263176");
    expect(str(S.binomialTest(7, 20, "1/2", { alternative: "less" }).pValue)).toBe("34495/262144");
    expect(str(S.binomialTest(7, 20, "1/2", { alternative: "greater" }).pValue)).toBe(str(R(1).subtract(S.binomial(20, "1/2").cdf(6))));
    // R: binom.test(3, 15, 0.5)$p.value = 0.03515625
    expect(str(S.binomialTest(3, 15, "0.5").pValue)).toBe("9/256");
    // an asymmetric null: P(X ≤ 2) + P(X = 10), the outcomes no more likely than 2
    expect(S.binomialTest(2, 10, "0.6").pValue.toNumber()).toBeCloseTo(0.0183411712, 12);
    const ci = t.interval();
    expect(ci.low.toNumber()).toBeCloseTo(0.1539092, 6);
    expect(ci.high.toNumber()).toBeCloseTo(0.5921885, 6);
    expect(() => S.binomialTest(5, 3)).toThrow(/more successes/);
    expect(() => S.binomialTest(1, 3, "1/2", { alternative: "sideways" })).toThrow(/alternative/);
  });
  test("intervals for a proportion", () => {
    const w = S.proportionInterval(7, 20, { method: "wilson" });
    expect(w.low.toNumber()).toBeCloseTo(0.1811919, 6);
    expect(w.high.toNumber()).toBeCloseTo(0.5671457, 6);
    const a = S.proportionInterval(7, 20, { method: "wald" });
    expect(a.low.toNumber()).toBeCloseTo(0.1409627, 6);
    expect(a.high.toNumber()).toBeCloseTo(0.5590373, 6);
    expect(str(S.proportionInterval(0, 10).low)).toBe("0");
    expect(S.proportionInterval(0, 10).high.toNumber()).toBeCloseTo(0.3084971, 6);
    expect(str(S.proportionInterval(10, 10).high)).toBe("1");
    expect(str(S.proportionInterval(0, 10, { method: "wilson" }).low)).toBe("0");
    expect(str(S.proportionInterval(10, 10, { method: "wilson" }).high)).toBe("1");
    // the defining property: at the lower end, seeing 7 or more has probability (1 − level)/2
    const low99 = S.proportionInterval(7, 20, { level: "0.99" }).low;
    expect(S.binomial(20, low99).sf(6).toNumber()).toBeCloseTo(0.005, 12);
    const high99 = S.proportionInterval(7, 20, { level: "0.99" }).high;
    expect(S.binomial(20, high99).cdf(7).toNumber()).toBeCloseTo(0.005, 12);
    expect(() => S.proportionInterval(7, 20, { method: "guess" })).toThrow(/unknown method/);
    expect(() => S.proportionInterval(7, 20, { level: 1 })).toThrow(/strictly between/);
    expect(() => S.proportionInterval(21, 20)).toThrow(/more successes/);
  });
  test("Fisher's test", () => {
    const t = S.fisherTest([[3, 1], [1, 3]]);
    expect(str(t.pValue)).toBe("17/35");
    expect(str(t.oddsRatio)).toBe("9");
    // R: fisher.test(matrix(c(8,1,2,5),2))$p.value = 0.03496503
    expect(str(S.fisherTest([[8, 2], [1, 5]]).pValue)).toBe("5/143");
    expect(str(S.fisherTest([[8, 2], [1, 5]], { alternative: "greater" }).pValue)).toBe("7/286");
    expect(S.fisherTest([[8, 2], [1, 5]], { alternative: "less" }).pValue.toNumber()).toBeCloseTo(0.9991259, 6);
    expect(S.fisherTest([[3, 0], [1, 3]]).oddsRatio).toBe(Infinity);
    expect(S.fisherTest([[0, 0], [1, 3]]).oddsRatio).toBe(null);
    expect(() => S.fisherTest([[1, 2, 3], [4, 5, 6]])).toThrow(/2×2/);
    expect(() => S.fisherTest([[0, 0], [0, 0]])).toThrow(/empty/);
  });
  test("permutation tests, exact and simulated", () => {
    const a = [12, 15, 11, 18, 14];
    const b = [9, 10, 13, 8, 11];
    const t = S.permutationTest(a, b);
    expect(t.kind).toBe("exact");
    expect(str(t.pValue)).toBe("1/21");
    expect(t.arrangements).toBe(252n);
    expect(str(t.statistic)).toBe("19/5");
    expect(str(S.permutationTest(a, b, { alternative: "greater" }).pValue)).toBe("1/42");
    expect(S.permutationTest(a, b, { alternative: "less" }).pValue.toNumber()).toBeGreaterThan(0.97);
    expect(S.permutationTest(a, b, { statistic: "median" }).kind).toBe("exact");
    expect(S.permutationTest(a, b, { statistic: (xs) => xs.reduce((m, v) => (v.compareTo(m) > 0 ? v : m)) }).kind).toBe("exact");
    const paired = S.permutationTest(a, b, { paired: true });
    expect(paired.arrangements).toBe(32n);
    expect(str(paired.pValue)).toBe("1/8");
    const sim = S.permutationTest(a, b, { limit: 10, resamples: 999, seed: 4 });
    expect(sim.kind).toBe("simulated");
    expect(sim.seed).toBe(4);
    expect(sim.pValueInterval.low.toNumber()).toBeLessThan(1 / 21);
    expect(sim.pValueInterval.high.toNumber()).toBeGreaterThan(1 / 21);
    expect(str(sim.pValue)).toBe(str(S.permutationTest(a, b, { limit: 10, resamples: 999, seed: 4 }).pValue));
    const simPaired = S.permutationTest(a, b, { paired: true, limit: 10, resamples: 499 });
    expect(simPaired.kind).toBe("simulated");
    expect(() => S.permutationTest(a, [1], { paired: true })).toThrow(/same length/);
    expect(() => S.permutationTest(a, b, { statistic: "mode" })).toThrow(/unknown statistic/);
    expect(() => S.permutationTest(["1:2"], [3])).toThrow(/exact data/);
  });
  test("a certified p-value can be too close to call", () => {
    const t = new S.TestResult({ test: "made up", kind: "certified", alternative: "two-sided", pValue: new RationalInterval(R("0.049"), R("0.051")) });
    expect(t.decision()).toBe("too close to call");
    expect(() => t.reject()).toThrow(/cannot be told apart/);
    expect(t.reject("0.1")).toBe(true);
  });
});
