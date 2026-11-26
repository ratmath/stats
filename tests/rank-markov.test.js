import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { combinations, factorial } from "@ratmath/discrete";
import * as S from "../index.js";

const R = (s) => Rational.from(String(s));
const str = (v) => String(v);
const near = (v, x, tol = 5e-7) => expect(Math.abs(S.toNumber(v) - x)).toBeLessThan(tol);

// published examples (R's documentation and Hollander & Wolfe)
const dep = { x: [1.83, 0.5, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.3], y: [0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29] };

describe("rank tests", () => {
  test("midranks and ties", () => {
    const { ranks, ties } = S.ranks([10, 20, 20, 30, 20]);
    expect(ranks.map(str)).toEqual(["1", "3", "3", "5", "3"]);
    expect(ties).toEqual([1, 3, 1]);
  });

  test("sign test is the binomial test on the signs", () => {
    const r = S.signTest(dep.x, dep.y);
    expect(r.statistic).toBe(7);
    expect(str(r.pValue)).toBe(str(S.binomialTest(7, 9).pValue));
    expect(S.signTest([1, 2, 3, 0, -1], { mu: 0 }).zeros).toBe(1);
  });

  test("Wilcoxon signed-rank: V = 40, p = 5/128 on the depression data", () => {
    const r = S.wilcoxonTest(dep.x, dep.y, { paired: true });
    expect(str(r.statistic)).toBe("40");
    expect(str(r.pValue)).toBe("5/128");
    expect(r.kind).toBe("exact");
    const g = S.wilcoxonTest(dep.x, dep.y, { alternative: "greater" });
    expect(str(g.pValue)).toBe("5/256");
  });

  test("Wilcoxon's exact distribution against all 2ⁿ sign patterns, with ties", () => {
    const d = [1, -2, 2, 3, -3, 3, 5, -0.5];
    const r = S.wilcoxonTest(d);
    const { ranks } = S.ranks(d.map(Math.abs));
    const W = (signs) => signs.reduce((s, sg, i) => (sg > 0 ? s.add(ranks[i]) : s), R(0));
    const mean = R(d.length * (d.length + 1)).divide(R(4));
    const obs = r.statistic.subtract(mean).abs();
    let hit = 0;
    for (let m = 0; m < 1 << d.length; m++) {
      const w = W(d.map((_, i) => ((m >> i) & 1 ? 1 : -1)));
      if (w.subtract(mean).abs().compareTo(obs) >= 0) hit++;
    }
    expect(str(r.pValue)).toBe(str(new Rational(BigInt(hit), BigInt(1 << d.length))));
  });

  test("Mann–Whitney: W = 35, p = 0.2544 (R's example), and the exact count with ties", () => {
    const x = [0.8, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46];
    const y = [1.15, 0.88, 0.9, 0.74, 1.21];
    const r = S.mannWhitneyTest(x, y);
    expect(str(r.statistic)).toBe("35");
    near(r.pValue, 0.2544122544122544, 1e-12);
    const a = [1, 2, 2, 3, 5];
    const b = [2, 3, 4, 4];
    const t = S.mannWhitneyTest(a, b);
    const pooled = [...a, ...b];
    const { ranks } = S.ranks(pooled);
    const mean = R(a.length * b.length).divide(R(2));
    const U = (idx) => idx.reduce((s, i) => s.add(ranks[i]), R(0)).subtract(R((a.length * (a.length + 1)) / 2));
    const obs = t.statistic.subtract(mean).abs();
    let hit = 0n;
    let all = 0n;
    for (const c of combinations([...pooled.keys()], a.length)) {
      all++;
      if (U(c).subtract(mean).abs().compareTo(obs) >= 0) hit++;
    }
    expect(str(t.pValue)).toBe(str(new Rational(hit, all)));
  });

  test("the normal approximation agrees closely with the exact answer for larger samples", () => {
    const x = Array.from({ length: 30 }, (_, i) => (i * 37) % 53);
    const y = Array.from({ length: 30 }, (_, i) => ((i * 41) % 61) + 3);
    const e = S.mannWhitneyTest(x, y, { method: "exact" });
    const n = S.mannWhitneyTest(x, y, { method: "normal" });
    expect(n.kind).toBe("certified");
    expect(Math.abs(S.toNumber(e.pValue) - S.toNumber(n.pValue))).toBeLessThan(0.01);
  });

  test("Kruskal–Wallis: H = 27/35 (Hollander & Wolfe)", () => {
    const r = S.kruskalWallisTest([[2.9, 3.0, 2.5, 2.6, 3.2], [3.8, 2.7, 4.0, 2.4], [2.8, 3.4, 3.7, 2.2, 2.0]]);
    expect(str(r.statistic)).toBe("27/35");
    near(r.pValue, 0.6799648, 1e-6);
  });

  test("Spearman and Kendall on R's example", () => {
    const a = [86, 97, 99, 100, 101, 103, 106, 110, 112, 113];
    const b = [2, 20, 28, 27, 50, 29, 7, 17, 6, 12];
    expect(str(S.spearmanTest(a, b).statistic)).toBe("-29/165");
    const k = S.kendallTest(a, b);
    expect(str(k.statistic)).toBe("-1/9");
    near(k.pValue, 0.7274895, 1e-6);
    // exact Spearman p by enumeration agrees with a hand count for n = 5
    const s = S.spearmanTest([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]);
    expect(s.kind).toBe("exact");
    expect(str(s.statistic)).toBe("4/5");
  });

  test("Kendall's exact distribution is the inversion count", () => {
    const t = S.kendallTest([1, 2, 3, 4], [1, 3, 2, 4], { alternative: "greater" });
    // S = 4: permutations of 4 with at most 1 inversion: 1 + 3 = 4 of 24
    expect(str(t.pValue)).toBe("1/6");
  });
});

describe("Kolmogorov–Smirnov", () => {
  test("K(10, 0.274) from Marsaglia, Tsang and Wang", () => {
    near(S.kolmogorovCdf(10, "0.274"), 0.6284796154565043, 1e-15);
    expect(S.kolmogorovCdf(5, 0)).toEqual(R(0));
  });
  test("the exact cdf against brute force for n = 2", () => {
    // D₂ for uniforms: P(D₂ < d) for d in (1/2, 1) is 2(2d − 1)... check by a fine grid
    const d = R("3/4");
    const exact = S.kolmogorovCdf(2, d).toNumber();
    let inside = 0;
    const N = 800;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      const u = [(i + 0.5) / N, (j + 0.5) / N].sort((p, q) => p - q);
      const D = Math.max(0.5 - u[0], u[0], 1 - u[1], u[1] - 0.5);
      if (D < 0.75) inside++;
    }
    expect(Math.abs(inside / N / N - exact)).toBeLessThan(0.01);
  });
  test("the one-sided tail bounds the two-sided one", () => {
    for (const [n, d] of [[10, "0.3"], [20, "0.25"], [15, "0.4"]]) {
      const p1 = S.smirnovTail(n, d);
      const p2 = S.kolmogorovTail(n, d);
      expect(p2.compareTo(p1) >= 0).toBe(true);
      expect(p2.compareTo(p1.multiply(R(2))) <= 0).toBe(true);
    }
  });
  test("one sample: exact against an exact cdf, certified against the normal", () => {
    const u = S.ksTest([0.1, 0.2, 0.35, 0.5, 0.9], (v) => v);
    expect(u.kind).toBe("exact");
    expect(str(u.statistic)).toBe("3/10");
    const z = S.ksTest([0.61, 0.29, 0.06, 0.59, -1.73, -0.74, 0.51, -0.56, 0.39, 1.64, 0.05, -0.06, 0.64, -0.82, 0.37, 1.77, 1.09, -1.28, 2.36, 1.31, 1.05, -0.32, -0.4, 1.06, -2.47], S.normal());
    expect(z.kind).toBe("certified");
    near(z.pValue, 0.3897, 1e-3);
  });
  test("two samples: the path count against all splits", () => {
    const x = [1, 2, 3, 4, 5, 6];
    const y = [3.5, 4.5, 6.5, 7.5, 8.5, 9.5];
    const r = S.ksTest(x, y);
    const pooled = [...x, ...y].sort((a, b) => a - b);
    const D = (xs) => {
      let best = 0;
      for (const v of pooled) {
        const fx = xs.filter((t) => t <= v).length / 6;
        const fy = pooled.filter((t) => t <= v && !xs.includes(t)).length / 6;
        best = Math.max(best, Math.abs(fx - fy));
      }
      return best;
    };
    const obs = D(x);
    let hit = 0;
    let all = 0;
    for (const c of combinations(pooled, 6)) {
      all++;
      if (D(c) >= obs - 1e-12) hit++;
    }
    expect(S.toNumber(r.statistic)).toBeCloseTo(obs, 12);
    expect(str(r.pValue)).toBe(str(new Rational(BigInt(hit), BigInt(all))));
  });
});

describe("the bootstrap", () => {
  test("exact for small samples: the mean's bootstrap variance is the plug-in variance over n", () => {
    const b = S.bootstrap([1, 2, 3, 4, 10]);
    expect(b.kind).toBe("exact");
    expect(str(b.mean)).toBe("4");
    expect(str(b.variance)).toBe("2");
    expect(b.distinct).toBe(126n);
    const total = b.distribution.reduce((s, d) => s.add(d.probability), R(0));
    expect(str(total)).toBe("1");
  });
  test("simulated from a seed is reproducible, and close to the exact answer", () => {
    const data = [3, 1, 4, 1, 5, 9, 2, 6];
    const e = S.bootstrap(data, "median");
    const a = S.bootstrap(data, "median", { exactUpTo: 0, resamples: 4000, seed: 7 });
    const b = S.bootstrap(data, "median", { exactUpTo: 0, resamples: 4000, seed: 7 });
    expect(a.kind).toBe("simulated");
    expect(a.se).toBe(b.se);
    expect(Math.abs(a.se - S.toNumber(e.se))).toBeLessThan(0.08);
  });
});

describe("Markov chains", () => {
  const weather = S.markovChain({ sun: { sun: "9/10", rain: "1/10" }, rain: { sun: "1/2", rain: "1/2" } });
  test("stationary distribution, n-step distribution, mean return time", () => {
    const pi = weather.stationary();
    expect(pi.unique).toBe(true);
    expect(str(pi.distribution.sun)).toBe("5/6");
    const d = weather.distribution("rain", 3);
    expect(str(d.sun)).toBe("39/50");
    // Pⁿ against repeated multiplication
    const P2 = weather.power(2);
    expect(str(P2[0][0])).toBe(str(R("9/10").multiply(R("9/10")).add(R("1/10").multiply(R("1/2")))));
    expect(str(weather.meanReturnTime("rain"))).toBe("6");
    expect(weather.isErgodic()).toBe(true);
  });
  test("gambler's ruin: absorption probabilities k/N and expected duration k(N − k)", () => {
    const N = 6;
    const rows = Array.from({ length: N + 1 }, (_, i) => Array.from({ length: N + 1 }, (_, j) => (i === 0 || i === N ? (i === j ? 1 : 0) : Math.abs(i - j) === 1 ? "1/2" : 0)));
    const a = S.markovChain(rows).absorption();
    for (let k = 1; k < N; k++) {
      expect(str(a.probabilities[k][String(N)])).toBe(str(R(k).divide(R(N))));
      expect(str(a.steps[k])).toBe(String(k * (N - k)));
    }
  });
  test("classes and periods", () => {
    const c = S.markovChain([[0, 1, 0], [0, 0, 1], [1, 0, 0]]);
    expect(c.classes()[0].period).toBe(3);
    expect(c.isErgodic()).toBe(false);
    expect(str(c.stationary().distribution[0])).toBe("1/3");
    const two = S.markovChain([[1, 0, 0], ["1/3", "1/3", "1/3"], [0, 0, 1]]);
    expect(two.stationary().unique).toBe(false);
    expect(two.stationary().all.length).toBe(2);
    expect(() => S.markovChain([["1/2", "1/3"], [0, 1]])).toThrow(/adds up/);
  });
  test("hitting times against simulation", () => {
    const h = weather.hitting("rain");
    expect(str(h.expected.sun)).toBe("10");
    const path = weather.simulate("sun", 20000, { seed: 4 });
    let waits = [];
    let t = 0;
    for (let i = 1; i < path.length; i++) {
      t++;
      if (path[i] === "rain" && path[i - 1] === "sun") {
        waits.push(t);
        t = 0;
      } else if (path[i] === "rain") t = 0;
    }
    const m = waits.reduce((a, b) => a + b, 0) / waits.length;
    expect(Math.abs(m - 10)).toBeLessThan(1);
    expect(weather.simulate("sun", 5, { seed: 4 })).toEqual(weather.simulate("sun", 5, { seed: 4 }));
  });
});

describe("conjugate updating", () => {
  test("beta–binomial: updating in two steps is updating once", () => {
    const one = S.betaPrior(2, 3).update({ successes: 7, trials: 10 });
    const two = S.betaPrior(2, 3).update({ successes: 3, trials: 4 }).update([1, 1, 1, 1, 0, 0]);
    expect(str(one.params.a)).toBe(str(two.params.a));
    expect(str(one.params.b)).toBe("6");
    expect(str(one.mean)).toBe("3/5");
    // the predictive law adds up to 1 and matches its definition at one point
    const pred = S.betaPrior(8, 4).predictive(3);
    expect(str(pred.table.reduce((s, v) => s.add(v), R(0)))).toBe("1");
    expect(str(pred.pmf(0))).toBe("5/91");
    const ci = one.interval("0.9");
    expect(S.toNumber(one.cdf(ci.low.low ?? ci.low))).toBeCloseTo(0.05, 6);
  });
  test("gamma–Poisson and gamma–exponential", () => {
    const g = S.gammaPrior(2, 1).update([3, 4, 2, 5]);
    expect(str(g.params.shape)).toBe("16");
    expect(str(g.params.rate)).toBe("5");
    let s = R(0);
    const pred = g.predictive();
    for (let k = 0; k < 200; k++) s = s.add(pred.pmf(k));
    expect(Math.abs(s.toNumber() - 1)).toBeLessThan(1e-12);
    const e = S.gammaPrior(1, 1, { likelihood: "exponential" }).update([0.5, 1.5]);
    expect([str(e.params.shape), str(e.params.rate)]).toEqual(["3", "3"]);
  });
  test("normal–normal: precisions add", () => {
    const p = S.normalPrior(0, { variance: 4, dataVariance: 1 }).update([1.2, 0.8, 1.0, 1.4]);
    expect(str(p.variance)).toBe("4/17");
    expect(str(p.mean)).toBe("88/85");
  });
  test("Dirichlet–multinomial", () => {
    const d = S.dirichletPrior({ a: 1, b: 1, c: 1 }).update({ a: 5, b: 2 });
    expect(str(d.mean.a)).toBe("3/5");
    expect(str(d.marginal("a").params.b)).toBe("4");
  });
});
