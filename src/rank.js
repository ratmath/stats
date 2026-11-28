// Rank tests, Kolmogorov–Smirnov and the bootstrap.
//
// A rank test uses only the order of the data, so its null distribution is a
// counting problem. Here it is counted: the p-value is an exact rational
// whenever the count is affordable (ties included, by midranks), and comes
// from the normal approximation, evaluated with certainty, for large samples.
// A result says which in `method`.
//
// Kolmogorov–Smirnov p-values are exact too: one sample by the matrix method
// of Marsaglia, Tsang and Wang, two samples by counting lattice paths. When
// the distribution's cdf is certified rather than exact, the statistic is an
// enclosure and so is the p-value. The bootstrap is simulated from a seed, or
// exact for very small samples, where every resample can be weighed.
import { Rational } from "@ratmath/core";
import { binomial as choose, factorial, rng as makeRng } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, TWO, HALF, Q, whole, probability, isQ, lowOf, highOf, interval, toNumber, show } from "./num.js";
import { values } from "./describe.js";
import { TestResult, binomialTest, proportionInterval } from "./tests.js";
import { chiSquare, fDist, Continuous } from "./continuous.js";
import { u01 } from "./continuous.js";

const ALTERNATIVES = new Set(["two-sided", "less", "greater"]);
const alternativeOf = (a) => {
  if (!ALTERNATIVES.has(a)) refuse("bad-argument", `the alternative must be "two-sided", "less" or "greater", got "${a}"`);
  return a;
};
const exactData = (data, what) => {
  const xs = values(data, what);
  if (xs.some((v) => !isQ(v))) refuse("needs-exact", `a rank test needs exact data, not intervals (${what})`);
  return xs;
};
const sign = (q) => (q.numerator > 0n ? 1 : q.numerator < 0n ? -1 : 0);
const clamp01 = (q) => (q.numerator < 0n ? ZERO : q.compareTo(ONE) > 0 ? ONE : q);
const METHODS = new Set(["auto", "exact", "normal"]);
const methodOf = (m) => {
  if (!METHODS.has(m)) refuse("bad-argument", `the method must be "auto", "exact" or "normal", got "${m}"`);
  return m;
};

/**
 * Midranks: the rank of each value in the list (1 for the smallest), ties
 * sharing the average of the ranks they span. Also the tie group sizes.
 */
export function ranks(data) {
  const xs = exactData(data, "the data");
  const order = xs.map((v, i) => i).sort((a, b) => xs[a].compareTo(xs[b]));
  const r = new Array(xs.length);
  const ties = [];
  for (let i = 0; i < order.length; ) {
    let j = i;
    while (j + 1 < order.length && xs[order[j + 1]].equals(xs[order[i]])) j++;
    const mid = Q(i + j + 2, 2);
    for (let k = i; k <= j; k++) r[order[k]] = mid;
    ties.push(j - i + 1);
    i = j + 1;
  }
  return { ranks: r, ties };
}
const tieTerm = (ties) => ties.reduce((s, t) => s + BigInt(t) ** 3n - BigInt(t), 0n);

/**
 * The p-value of a statistic S from its exact null distribution, given as a
 * map from twice the statistic (an integer) to a count, out of `total`.
 * Two-sided: P(|S − E S| ≥ |s − E S|).
 */
function tailFromCounts(counts, total, s2, alt) {
  let mean2 = ZERO;
  for (const [v, c] of counts) mean2 = mean2.add(new Rational(v * c));
  mean2 = mean2.divide(new Rational(total));
  const far = (v) => new Rational(v).subtract(mean2).abs();
  const d = far(s2);
  let hit = 0n;
  for (const [v, c] of counts) {
    if (alt === "greater" ? v >= s2 : alt === "less" ? v <= s2 : far(v).compareTo(d) >= 0) hit += c;
  }
  return new Rational(hit, total);
}
/** The normal approximation: z² = (S − E S ∓ correction)²/Var S, exact; the p-value certified. */
function normalP(S, mean, variance, alt, { correct = true, step = ONE, bits = 64 } = {}) {
  let d = S.subtract(mean);
  const s = sign(d);
  if (correct && s !== 0) {
    const c = step.divide(TWO);
    d = d.abs().compareTo(c) <= 0 ? ZERO : s > 0 ? d.subtract(c) : d.add(c);
  }
  if (variance.numerator === 0n) refuse("no-variation", "every value is tied, so the ranks carry no information");
  const z2 = d.multiply(d).divide(variance);
  const two = chiSquare(1).sf(z2, { bits });
  let pValue;
  if (alt === "two-sided") pValue = two;
  else {
    const towards = alt === "greater" ? s > 0 : s < 0;
    const half = isQ(two) ? two.divide(TWO) : interval(two.low.divide(TWO), two.high.divide(TWO));
    pValue = s === 0 ? HALF : towards ? half : isQ(half) ? ONE.subtract(half) : interval(ONE.subtract(half.high), ONE.subtract(half.low));
  }
  return { z2, pValue };
}

/**
 * The sign test: under H₀ the median of x (or of x − y, paired) is mu, so
 * each non-zero difference is positive with probability 1/2. Zeros are
 * dropped. Exact, by the binomial distribution.
 */
export function signTest(x, y, { mu = 0, alternative = "two-sided" } = {}) {
  if (y !== undefined && !Array.isArray(y) && typeof y === "object") return signTest(x, undefined, y);
  const d = differences(x, y, mu);
  const nz = d.filter((v) => v.numerator !== 0n);
  if (!nz.length) refuse("no-data", "every difference is zero");
  const plus = nz.filter((v) => v.numerator > 0n).length;
  const b = binomialTest(plus, nz.length, HALF, { alternative: alternativeOf(alternative) });
  return new TestResult({
    test: y === undefined ? "sign test" : "paired sign test", kind: "exact", method: "exact", statistic: plus, n: nz.length,
    zeros: d.length - nz.length, alternative: b.alternative, pValue: b.pValue, null: Q(toQ(mu, "mu").numerator, toQ(mu, "mu").denominator),
  });
}
function differences(x, y, mu) {
  const m = toQ(mu, "mu");
  const xs = exactData(x, "x");
  if (y === undefined) return xs.map((v) => v.subtract(m));
  const ys = exactData(y, "y");
  if (ys.length !== xs.length) refuse("bad-argument", "paired data need lists of the same length");
  return xs.map((v, i) => v.subtract(ys[i]).subtract(m));
}

/**
 * Wilcoxon's signed-rank test of H₀: the differences (x − mu, or x − y − mu
 * when paired) are symmetric about 0. W is the sum of the ranks of |d| over
 * the positive d; zeros are dropped and ties take midranks.
 *   method "exact"   counts the 2ⁿ sign patterns (by dynamic programming)
 *   method "normal"  z from the tie-corrected variance, with continuity correction
 *   method "auto"    exact up to n = 60
 * Also the Hodges–Lehmann estimate: the median of the Walsh averages.
 */
export function wilcoxonTest(x, y, options = {}) {
  if (y !== undefined && !Array.isArray(y) && typeof y === "object" && !(Symbol.iterator in y)) return wilcoxonTest(x, undefined, y);
  const { mu = 0, alternative = "two-sided", method = "auto", correct = true, bits = 64 } = options;
  const alt = alternativeOf(alternative);
  const how = methodOf(method);
  const d = differences(x, y, mu).filter((v) => v.numerator !== 0n);
  const n = d.length;
  if (!n) refuse("no-data", "every difference is zero");
  const { ranks: r, ties } = ranks(d.map((v) => v.abs()));
  let W = ZERO;
  d.forEach((v, i) => {
    if (v.numerator > 0n) W = W.add(r[i]);
  });
  const mean = Q(n * (n + 1), 4);
  const variance = Q(BigInt(n) * BigInt(n + 1) * BigInt(2 * n + 1), 24n).subtract(new Rational(tieTerm(ties), 48n));
  const walsh = [];
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) walsh.push(d[i].add(d[j]).divide(TWO));
  const estimate = medianOf(walsh).add(toQ(mu, "mu"));
  const base = { test: y === undefined ? "Wilcoxon signed-rank test" : "Wilcoxon signed-rank test (paired)", statistic: W, n, alternative: alt, mean, variance, estimate };
  if (how === "exact" || (how === "auto" && n <= 60)) {
    // twice each midrank is a whole number; count sign patterns by the sum of 2·rank over the positives
    let counts = new Map([[0n, 1n]]);
    for (const rk of r) {
      const w = rk.numerator * (2n / rk.denominator);
      const next = new Map(counts);
      for (const [s, c] of counts) next.set(s + w, (next.get(s + w) ?? 0n) + c);
      counts = next;
    }
    const pValue = tailFromCounts(counts, 1n << BigInt(n), W.numerator * (2n / W.denominator), alt);
    return new TestResult({ ...base, kind: "exact", method: "exact", pValue });
  }
  const { z2, pValue } = normalP(W, mean, variance, alt, { correct, bits });
  return new TestResult({ ...base, kind: "certified", method: "normal", z2, pValue });
}
const medianOf = (list) => {
  const s = [...list].sort((a, b) => a.compareTo(b));
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : s[m - 1].add(s[m]).divide(TWO);
};

/**
 * The Mann–Whitney (Wilcoxon rank-sum) test of H₀: x and y come from the same
 * distribution. U counts the pairs with x above y (ties count 1/2). "greater"
 * means x tends to be larger.
 *   method "exact"   counts every way to choose which pooled values are x's,
 *                    with the observed midranks (so ties are exact too)
 *   method "normal"  tie-corrected variance, continuity correction
 *   method "auto"    exact when both groups are at most 50 and together 80
 * Also the Hodges–Lehmann shift estimate: the median of the differences x − y.
 */
export function mannWhitneyTest(x, y, { alternative = "two-sided", method = "auto", correct = true, bits = 64 } = {}) {
  const alt = alternativeOf(alternative);
  const how = methodOf(method);
  const xs = exactData(x, "x");
  const ys = exactData(y, "y");
  const m = xs.length;
  const n = ys.length;
  const N = m + n;
  const { ranks: r, ties } = ranks([...xs, ...ys]);
  const R1 = r.slice(0, m).reduce((s, v) => s.add(v), ZERO);
  const U = R1.subtract(Q(m * (m + 1), 2));
  const mean = Q(m * n, 2);
  const variance = Q(m * n, 12).multiply(Q(N + 1).subtract(new Rational(tieTerm(ties), BigInt(N) * BigInt(N - 1) || 1n)));
  const diffs = [];
  for (const a of xs) for (const b of ys) diffs.push(a.subtract(b));
  const base = { test: "Mann–Whitney test", statistic: U, rankSum: R1, m, n, alternative: alt, mean, variance, estimate: medianOf(diffs) };
  if (how === "exact" || (how === "auto" && m <= 50 && n <= 50 && N <= 80)) {
    // counts[k]: for k items chosen, a map from the sum of their doubled ranks to the number of ways
    let counts = Array.from({ length: m + 1 }, () => new Map());
    counts[0].set(0n, 1n);
    for (const rk of r) {
      const w = rk.numerator * (2n / rk.denominator);
      for (let k = m; k >= 1; k--) {
        const from = counts[k - 1];
        if (!from.size) continue;
        const to = counts[k];
        for (const [s, c] of from) to.set(s + w, (to.get(s + w) ?? 0n) + c);
      }
    }
    const shift = BigInt(m * (m + 1)); // 2·m(m+1)/2
    const uCounts = new Map([...counts[m]].map(([s, c]) => [s - shift, c]));
    const pValue = tailFromCounts(uCounts, choose(N, m), U.numerator * (2n / U.denominator), alt);
    return new TestResult({ ...base, kind: "exact", method: "exact", pValue });
  }
  const { z2, pValue } = normalP(U, mean, variance, alt, { correct, bits });
  return new TestResult({ ...base, kind: "certified", method: "normal", z2, pValue });
}

/**
 * The Kruskal–Wallis test of H₀: k groups come from the same distribution.
 * H is exact (tie-corrected); the p-value is from chi-square with k − 1
 * degrees of freedom, certified.
 */
export function kruskalWallisTest(groups, { bits = 64 } = {}) {
  if (!Array.isArray(groups) || groups.length < 2) refuse("bad-argument", "kruskalWallisTest needs a list of two or more groups");
  const gs = groups.map((g, i) => exactData(g, `group ${i + 1}`));
  const all = gs.flat();
  const N = all.length;
  const { ranks: r, ties } = ranks(all);
  let at = 0;
  let S = ZERO;
  const meanRanks = gs.map((g) => {
    const R = r.slice(at, at + g.length).reduce((s, v) => s.add(v), ZERO);
    at += g.length;
    S = S.add(R.multiply(R).divide(Q(g.length)));
    return R.divide(Q(g.length));
  });
  let H = Q(12).divide(Q(N * (N + 1))).multiply(S).subtract(Q(3 * (N + 1)));
  const C = ONE.subtract(new Rational(tieTerm(ties), BigInt(N) ** 3n - BigInt(N)));
  if (C.numerator === 0n) refuse("no-variation", "every value is tied");
  H = H.divide(C);
  const df = gs.length - 1;
  return new TestResult({
    test: "Kruskal–Wallis test", kind: "certified", method: "chi-square", statistic: H, df, meanRanks, alternative: "two-sided",
    pValue: chiSquare(df).sf(H, { bits }),
  });
}

/**
 * Spearman's rank correlation: Pearson's correlation of the midranks. ρ² is
 * exact; ρ is exact when there are no ties and certified otherwise.
 *   method "exact"   all n! orderings of y against x (n ≤ 9)
 *   method "normal"  t² = ρ²(n − 2)/(1 − ρ²) against F(1, n − 2), certified
 */
export function spearmanTest(x, y, { alternative = "two-sided", method = "auto", bits = 64 } = {}) {
  const alt = alternativeOf(alternative);
  const how = methodOf(method);
  const xs = exactData(x, "x");
  const ys = exactData(y, "y");
  if (xs.length !== ys.length) refuse("bad-argument", "x and y need the same length");
  const n = xs.length;
  if (n < 3) refuse("no-data", "a rank correlation needs at least three pairs");
  const rx = ranks(xs).ranks;
  const ry = ranks(ys).ranks;
  const c = (Rn) => {
    const m = Rn.reduce((s, v) => s.add(v), ZERO).divide(Q(n));
    return Rn.map((v) => v.subtract(m));
  };
  const cx = c(rx);
  const cy = c(ry);
  const dot = (a, b) => a.reduce((s, v, i) => s.add(v.multiply(b[i])), ZERO);
  const sxx = dot(cx, cx);
  const syy = dot(cy, cy);
  if (sxx.numerator === 0n || syy.numerator === 0n) refuse("no-variation", "one of the variables is constant");
  const sxy = dot(cx, cy);
  const rho2 = sxy.multiply(sxy).divide(sxx.multiply(syy));
  const s = sign(sxy);
  const tiesFree = sxx.equals(syy) && sxx.equals(Q(BigInt(n) * BigInt(n * n - 1), 12n));
  let rho;
  if (sxx.equals(syy)) rho = sxy.divide(sxx);
  else {
    const root = sqrtQ(rho2, bits);
    rho = s < 0 ? interval(root.high.negate(), root.low.negate()) : root;
  }
  const base = { test: "Spearman's rank correlation test", statistic: rho, rho2, n, alternative: alt };
  if (how === "exact" || (how === "auto" && n <= 9)) {
    if (n > 10) refuse("too-large", "the exact Spearman test counts n! orderings; use method \"normal\" for n above 10");
    const perms = permutationsOf(n);
    let hit = 0n;
    let total = 0n;
    for (const p of perms) {
      const v = dot(cx, p.map((i) => cy[i]));
      total++;
      const ok = alt === "greater" ? v.compareTo(sxy) >= 0 : alt === "less" ? v.compareTo(sxy) <= 0 : v.abs().compareTo(sxy.abs()) >= 0;
      if (ok) hit++;
    }
    return new TestResult({ ...base, kind: "exact", method: "exact", pValue: new Rational(hit, total), tiesFree });
  }
  if (rho2.equals(ONE)) return new TestResult({ ...base, kind: "exact", method: "normal", pValue: alt === "two-sided" || (alt === "greater") === s > 0 ? ZERO : ONE });
  const t2 = rho2.multiply(Q(n - 2)).divide(ONE.subtract(rho2));
  const two = fDist(1, n - 2).sf(t2, { bits });
  return new TestResult({ ...base, kind: "certified", method: "t", t2, pValue: oneSide(two, s, alt) });
}
const oneSide = (two, s, alt) => {
  if (alt === "two-sided") return two;
  if (s === 0) return HALF;
  const towards = alt === "greater" ? s > 0 : s < 0;
  const half = isQ(two) ? two.divide(TWO) : interval(two.low.divide(TWO), two.high.divide(TWO));
  return towards ? half : isQ(half) ? ONE.subtract(half) : interval(ONE.subtract(half.high), ONE.subtract(half.low));
};
function* permutationsOf(n) {
  const a = Array.from({ length: n }, (_, i) => i);
  const c = new Array(n).fill(0);
  yield a.slice();
  for (let i = 0; i < n; ) {
    if (c[i] < i) {
      const j = i % 2 ? c[i] : 0;
      [a[j], a[i]] = [a[i], a[j]];
      yield a.slice();
      c[i]++;
      i = 0;
    } else c[i++] = 0;
  }
}
/** An enclosure of √q for a rational q ≥ 0, about `bits` bits wide. */
function sqrtQ(q, bits) {
  const scale = 1n << BigInt(bits);
  const target = (q.numerator * scale * scale) / q.denominator;
  let r = sqrtBig(target);
  return interval(new Rational(r, scale), new Rational(r + 1n, scale));
}
function sqrtBig(n) {
  if (n < 2n) return n;
  let x = BigInt(Math.floor(Math.sqrt(Number(n))));
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x && x * x <= n && (x + 1n) * (x + 1n) > n) return x;
    x = y;
  }
}

/**
 * Kendall's τ: S = concordant − discordant pairs, and τ_b = S/√(n₀ − n₁)(n₀ − n₂)
 * with the tie corrections (τ_a = τ_b without ties).
 *   method "exact"   the distribution of S with no ties, by counting
 *                    permutations by inversions (n ≤ 60)
 *   method "normal"  z² = S²/Var S with the tie-corrected variance, certified
 */
export function kendallTest(x, y, { alternative = "two-sided", method = "auto", bits = 64 } = {}) {
  const alt = alternativeOf(alternative);
  const how = methodOf(method);
  const xs = exactData(x, "x");
  const ys = exactData(y, "y");
  if (xs.length !== ys.length) refuse("bad-argument", "x and y need the same length");
  const n = xs.length;
  if (n < 3) refuse("no-data", "a rank correlation needs at least three pairs");
  let S = 0;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) S += xs[i].compareTo(xs[j]) * ys[i].compareTo(ys[j]);
  const tx = ranks(xs).ties;
  const ty = ranks(ys).ties;
  const n0 = (n * (n - 1)) / 2;
  const pairs = (t) => t.reduce((s, k) => s + (k * (k - 1)) / 2, 0);
  const n1 = pairs(tx);
  const n2 = pairs(ty);
  const denom2 = Q((n0 - n1) * (n0 - n2));
  if (denom2.numerator === 0n) refuse("no-variation", "one of the variables is constant");
  const noTies = n1 === 0 && n2 === 0;
  const tau = noTies ? Q(S, n0) : (() => {
    const r = sqrtQ(Q(S * S).divide(denom2), bits);
    return S < 0 ? interval(r.high.negate(), r.low.negate()) : r;
  })();
  const base = { test: "Kendall's rank correlation test", statistic: tau, S, n, alternative: alt };
  if (noTies && (how === "exact" || (how === "auto" && n <= 60))) {
    // the number of permutations of n with k inversions; S = n0 − 2·inversions
    let c = [1n];
    for (let m = 2; m <= n; m++) {
      const next = new Array(c.length + m - 1).fill(0n);
      for (let k = 0; k < c.length; k++) if (c[k]) for (let j = 0; j < m; j++) next[k + j] += c[k];
      c = next;
    }
    const counts = new Map(c.map((cnt, inv) => [BigInt(2 * (n0 - 2 * inv)), cnt]));
    return new TestResult({ ...base, kind: "exact", method: "exact", pValue: tailFromCounts(counts, factorial(n), BigInt(2 * S), alt) });
  }
  const v = (t, f) => t.reduce((s, k) => s + f(k), 0);
  const big = (k) => k * (k - 1) * (2 * k + 5);
  let variance = Q(big(n) - v(tx, big) - v(ty, big), 18);
  variance = variance.add(Q(v(tx, (k) => k * (k - 1) * (k - 2)) * v(ty, (k) => k * (k - 1) * (k - 2)), 9 * n * (n - 1) * (n - 2)));
  variance = variance.add(Q(v(tx, (k) => k * (k - 1)) * v(ty, (k) => k * (k - 1)), 2 * n * (n - 1)));
  const { z2, pValue } = normalP(Q(S), ZERO, variance, alt, { correct: false, bits });
  return new TestResult({ ...base, kind: "certified", method: "normal", z2, variance, pValue });
}

// ---- Kolmogorov–Smirnov ------------------------------------------------------------------------

/**
 * P(Dₙ < d) for the one-sample Kolmogorov statistic, exactly, by the method
 * of Marsaglia, Tsang and Wang (2003): n!/nⁿ times an entry of Hⁿ, with H a
 * (2k − 1)-square matrix of rationals, k = ⌈nd⌉.
 */
export function kolmogorovCdf(n, d) {
  const N = whole(n, "n", { min: 1 });
  const D = toQ(d, "d");
  if (D.numerator <= 0n) return ZERO;
  if (D.compareTo(ONE) >= 0) return ONE;
  const nd = D.multiply(Q(N));
  // k = ⌈nd⌉
  const K = Math.max(1, Number(nd.denominator === 1n ? nd.numerator : BigInt(nd.floor()) + 1n));
  const h = Q(K).subtract(nd);
  const m = 2 * K - 1;
  const fact = [ONE];
  for (let i = 1; i <= m + 1; i++) fact.push(fact[i - 1].multiply(Q(i)));
  const H = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i - j + 1 >= 0 ? ONE : ZERO)));
  for (let i = 0; i < m; i++) {
    H[i][0] = H[i][0].subtract(h.pow(BigInt(i + 1)));
    H[m - 1][i] = H[m - 1][i].subtract(h.pow(BigInt(m - i)));
  }
  const t = TWO.multiply(h).subtract(ONE);
  if (t.numerator > 0n) H[m - 1][0] = H[m - 1][0].add(t.pow(BigInt(m)));
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) if (i - j + 1 > 0) H[i][j] = H[i][j].divide(fact[i - j + 1]);
  const mul = (A, B) => A.map((row) => B[0].map((_, j) => {
    let s = ZERO;
    for (let t2 = 0; t2 < m; t2++) if (row[t2].numerator !== 0n && B[t2][j].numerator !== 0n) s = s.add(row[t2].multiply(B[t2][j]));
    return s;
  }));
  let R = null;
  let P = H;
  let e = N;
  while (e > 0) {
    if (e & 1) R = R ? mul(R, P) : P;
    e >>= 1;
    if (e) P = mul(P, P);
  }
  const v = R[K - 1][K - 1].multiply(new Rational(factorial(N), BigInt(N) ** BigInt(N)));
  return clamp01(v);
}

/**
 * P(Dₙ⁺ ≥ d), the one-sided tail, exactly (Smirnov; Birnbaum and Tingey):
 *   d · Σ_{j=0}^{⌊n(1−d)⌋} C(n, j) (1 − d − j/n)^{n−j} (d + j/n)^{j−1}.
 */
export function smirnovTail(n, d) {
  const N = whole(n, "n", { min: 1 });
  const D = toQ(d, "d");
  if (D.numerator <= 0n) return ONE;
  if (D.compareTo(ONE) >= 0) return D.equals(ONE) ? ONE.divide(Q(N) .pow(BigInt(N))) : ZERO;
  const J = Number(ONE.subtract(D).multiply(Q(N)).floor());
  let s = ZERO;
  for (let j = 0; j <= J; j++) {
    const a = ONE.subtract(D).subtract(Q(j, N));
    const b = D.add(Q(j, N));
    s = s.add(new Rational(choose(N, j)).multiply(a.pow(BigInt(N - j))).multiply(j === 0 ? ONE.divide(b) : b.pow(BigInt(j - 1))));
  }
  return clamp01(D.multiply(s));
}
/**
 * The two-sided tail P(Dₙ ≥ d). Dₙ ≥ d exactly when Dₙ⁺ ≥ d or Dₙ⁻ ≥ d,
 * which have the same law, so P⁺ ≤ P(Dₙ ≥ d) ≤ 2P⁺. Far in the tail
 * (2P⁺ ≤ 1/1000) that enclosure is returned, since it decides any usual
 * level; otherwise the exact value by the matrix method.
 */
export function kolmogorovTail(n, d) {
  const N = whole(n, "n", { min: 1 });
  const D = toQ(d, "d");
  const p = smirnovTail(N, D);
  const hi = TWO.multiply(p);
  if (hi.compareTo(Q(1, 1000)) > 0) return ONE.subtract(kolmogorovCdf(N, D));
  return p.numerator === 0n ? ZERO : interval(p, hi);
}

/**
 * The Kolmogorov–Smirnov test.
 *   ksTest(x, dist)   one sample against a distribution: a Continuous (its
 *                     certified cdf) or a function x ↦ F(x) giving exact
 *                     rationals. Two-sided only. The p-value is exact when F
 *                     is, otherwise an enclosure; for n above `exactUpTo`
 *                     (default 100) it is simulated from `seed`.
 *   ksTest(x, y)      two samples: the largest gap between the empirical
 *                     cdfs, with the exact p-value of the permutation
 *                     distribution (ties respected), by counting lattice paths.
 *                     "greater" means x tends to be larger (its cdf lies below).
 */
export function ksTest(x, other, { alternative = "two-sided", exactUpTo = 100, resamples = 2000, seed = 1, bits = 96 } = {}) {
  const alt = alternativeOf(alternative);
  const xs = exactData(x, "x").sort((a, b) => a.compareTo(b));
  const n = xs.length;
  if (Array.isArray(other) || (other && typeof other[Symbol.iterator] === "function" && !(other instanceof Continuous))) return ksTwo(xs, exactData(other, "y").sort((a, b) => a.compareTo(b)), alt);
  if (alt !== "two-sided") refuse("bad-argument", "the one-sample test here is two-sided");
  const F = other instanceof Continuous ? (v) => other.cdf(v, { bits }) : typeof other === "function" ? (v) => toQ(other(v), "F(x)") : refuse("bad-argument", "ksTest needs a second sample, a Continuous distribution, or a cdf function");
  let dLo = ZERO;
  let dHi = ZERO;
  xs.forEach((v, i) => {
    const f = F(v);
    const up = Q(i + 1, n);
    const down = Q(i, n);
    // the gaps i/n − F and F − (i−1)/n, as intervals
    const a = [up.subtract(highOf(f)), up.subtract(lowOf(f))];
    const b = [lowOf(f).subtract(down), highOf(f).subtract(down)];
    for (const [lo, hi] of [a, b]) {
      if (lo.compareTo(dLo) > 0) dLo = lo;
      if (hi.compareTo(dHi) > 0) dHi = hi;
    }
  });
  const D = dLo.equals(dHi) ? dLo : interval(dLo, dHi);
  const base = { test: "Kolmogorov–Smirnov test", statistic: D, n, alternative: alt, against: other instanceof Continuous ? String(other.name ?? "the distribution") : "F" };
  if (n <= exactUpTo) {
    const pAt = (d) => kolmogorovTail(n, d);
    if (isQ(D)) {
      const p = pAt(D);
      return new TestResult({ ...base, kind: isQ(p) ? "exact" : "certified", method: isQ(p) ? "exact" : "one-sided bound", pValue: p });
    }
    // P(Dₙ ≥ d) falls as d grows. The ends of D are rounded outward to a
    // multiple of 1/(n·2²⁰) first: the enclosure stays valid and the matrices stay small.
    const S = BigInt(n) << 20n;
    const down = (q) => new Rational((q.numerator * S) / q.denominator, S);
    const up = (q) => new Rational(-((-q.numerator * S) / q.denominator) + ((q.numerator * S) % q.denominator === 0n ? 0n : q.numerator > 0n ? 1n : 0n), S);
    const lo = lowOf(pAt(up(dHi)));
    const hi = highOf(pAt(down(dLo)));
    return new TestResult({ ...base, kind: "certified", method: "exact", pValue: lo.equals(hi) ? lo : interval(lo, hi) });
  }
  // simulated: Dₙ has the same law for every continuous F, so draw uniforms
  const B = whole(resamples, "resamples", { min: 1 });
  const g = makeRng(seed).fork("Kolmogorov–Smirnov");
  const dNum = toNumber(dLo);
  let hits = 0;
  for (let r = 0; r < B; r++) {
    const u = Array.from({ length: n }, () => u01(g)).sort((a, b) => a - b);
    let dd = 0;
    u.forEach((v, i) => (dd = Math.max(dd, (i + 1) / n - v, v - i / n)));
    if (dd >= dNum) hits++;
  }
  const ci = proportionInterval(hits, B, { level: "0.99" });
  return new TestResult({ ...base, kind: "simulated", method: "simulated", seed, resamples: B, hits, pValue: Q(hits + 1, B + 1), pValueInterval: { low: ci.low, high: ci.high, level: ci.level } });
}
function ksTwo(xs, ys, alt) {
  const m = xs.length;
  const n = ys.length;
  // walk the pooled order; check the gap only where the value changes
  const pooled = [...xs.map((v) => [v, 0]), ...ys.map((v) => [v, 1])].sort((a, b) => a[0].compareTo(b[0]));
  const gap = (i, j) => Q(i, m).subtract(Q(j, n)); // F_x − F_y after i x's and j y's
  let i = 0;
  let j = 0;
  let Dp = ZERO; // max of F_y − F_x: x larger
  let Dm = ZERO;
  const checkpoints = new Set(); // positions (i + j) after which the value changes
  for (let k = 0; k < pooled.length; k++) {
    if (pooled[k][1] === 0) i++;
    else j++;
    if (k + 1 === pooled.length || !pooled[k + 1][0].equals(pooled[k][0])) {
      checkpoints.add(k + 1);
      const g = gap(i, j);
      if (g.negate().compareTo(Dp) > 0) Dp = g.negate();
      if (g.compareTo(Dm) > 0) Dm = g;
    }
  }
  const D = alt === "greater" ? Dp : alt === "less" ? Dm : Dp.compareTo(Dm) >= 0 ? Dp : Dm;
  // count paths from (0,0) to (m,n) that stay strictly inside the band at every checkpoint
  const inside = (a, b) => {
    const g = gap(a, b);
    const v = alt === "greater" ? g.negate() : alt === "less" ? g : g.abs();
    return v.compareTo(D) < 0;
  };
  let row = new Array(n + 1).fill(0n);
  for (let a = 0; a <= m; a++) {
    const next = new Array(n + 1).fill(0n);
    for (let b = 0; b <= n; b++) {
      let c = a === 0 && b === 0 ? 1n : (a > 0 ? row[b] : 0n) + (b > 0 ? next[b - 1] : 0n);
      if (c && a + b > 0 && checkpoints.has(a + b) && !inside(a, b)) c = 0n;
      next[b] = c;
    }
    row = next;
  }
  const total = choose(m + n, m);
  return new TestResult({
    test: "two-sample Kolmogorov–Smirnov test", kind: "exact", method: "exact", statistic: D, m, n, alternative: alt,
    pValue: ONE.subtract(new Rational(row[n], total)),
  });
}

// ---- the bootstrap -----------------------------------------------------------------------------

const STATS = {
  mean: (xs) => xs.reduce((s, v) => s.add(v), ZERO).divide(Q(xs.length)),
  median: medianOf,
  variance(xs) {
    const m = STATS.mean(xs);
    return xs.length < 2 ? ZERO : xs.reduce((s, v) => s.add(v.subtract(m).multiply(v.subtract(m))), ZERO).divide(Q(xs.length - 1));
  },
};
/**
 * The bootstrap: the sampling distribution of a statistic, estimated by
 * resampling the data with replacement.
 *   bootstrap(data, statistic, { resamples, seed, level, method, exactUpTo })
 * statistic: "mean", "median", "variance", or a function of a list of
 * Rationals returning a number or Rational.
 *   method "percentile" (the default) or "basic" for the interval.
 * When the number of distinct resamples C(2n − 1, n) is at most `exactUpTo`
 * (default 20,000), every resample is weighed by its multinomial probability
 * and the answer is exact: the bootstrap distribution itself, not a
 * simulation of it. Otherwise `resamples` draws from `seed`.
 * Returns { estimate, mean, bias, se, interval, kind, ... }.
 */
export function bootstrap(data, statistic = "mean", { resamples = 2000, seed = 1, level = "0.95", method = "percentile", exactUpTo = 20000, bits = 64 } = {}) {
  const xs = exactData(data, "the data");
  const stat = typeof statistic === "function" ? (l) => toQ(statistic(l), "the statistic") : STATS[statistic];
  if (!stat) refuse("bad-argument", `unknown statistic "${statistic}"; use "mean", "median", "variance" or a function`);
  if (method !== "percentile" && method !== "basic") refuse("bad-argument", `the method must be "percentile" or "basic"`);
  const conf = probability(level, "the level");
  const alpha = ONE.subtract(conf).divide(TWO);
  const n = xs.length;
  const estimate = stat(xs);
  const distinct = choose(2 * n - 1, n);
  const toInterval = (lo, hi) => (method === "basic" ? { low: TWO.multiply(estimate).subtract(hi), high: TWO.multiply(estimate).subtract(lo) } : { low: lo, high: hi });
  if (distinct <= BigInt(exactUpTo)) {
    // every multiset of size n from the n positions, with probability n!/(∏cᵢ!·nⁿ)
    const dist = [];
    const counts = new Array(n).fill(0);
    const denom = BigInt(n) ** BigInt(n);
    const nf = factorial(n);
    const go = (i, left) => {
      if (i === n - 1) {
        counts[i] = left;
        const sample = [];
        let w = nf;
        counts.forEach((c, k) => {
          for (let t = 0; t < c; t++) sample.push(xs[k]);
          w /= factorial(c);
        });
        dist.push([stat(sample), new Rational(w, denom)]);
        return;
      }
      for (let c = left; c >= 0; c--) {
        counts[i] = c;
        go(i + 1, left - c);
      }
    };
    go(0, n);
    dist.sort((a, b) => a[0].compareTo(b[0]));
    // one entry per value
    for (let i = dist.length - 1; i > 0; i--) {
      if (dist[i][0].equals(dist[i - 1][0])) {
        dist[i - 1][1] = dist[i - 1][1].add(dist[i][1]);
        dist.splice(i, 1);
      }
    }
    let mean = ZERO;
    for (const [v, p] of dist) mean = mean.add(v.multiply(p));
    let variance = ZERO;
    for (const [v, p] of dist) variance = variance.add(v.subtract(mean).multiply(v.subtract(mean)).multiply(p));
    const quantileAt = (q) => {
      let s = ZERO;
      for (const [v, p] of dist) {
        s = s.add(p);
        if (s.compareTo(q) >= 0) return v;
      }
      return dist[dist.length - 1][0];
    };
    const lo = quantileAt(alpha);
    const hi = quantileAt(ONE.subtract(alpha));
    return new Bootstrap({
      kind: "exact", method, statistic: typeof statistic === "string" ? statistic : "custom", n, estimate, mean, bias: mean.subtract(estimate), variance,
      se: sqrtQ(variance, bits), interval: { ...toInterval(lo, hi), level: conf }, distinct, distribution: dist.map(([value, probability]) => ({ value, probability })),
    });
  }
  const B = whole(resamples, "resamples", { min: 2 });
  const g = makeRng(seed).fork("bootstrap");
  const reps = new Array(B);
  for (let r = 0; r < B; r++) {
    const sample = new Array(n);
    for (let i = 0; i < n; i++) sample[i] = xs[g.int(n)];
    reps[r] = stat(sample);
  }
  reps.sort((a, b) => a.compareTo(b));
  const nums = reps.map((v) => v.toNumber());
  const meanN = nums.reduce((a, b) => a + b, 0) / B;
  const sdN = Math.sqrt(nums.reduce((s, v) => s + (v - meanN) ** 2, 0) / (B - 1));
  const at = (q) => reps[Math.min(B - 1, Math.max(0, Math.ceil(q.multiply(Q(B)).toNumber()) - 1))];
  const lo = at(alpha);
  const hi = at(ONE.subtract(alpha));
  return new Bootstrap({
    kind: "simulated", method, statistic: typeof statistic === "string" ? statistic : "custom", n, estimate, seed, resamples: B,
    mean: meanN, bias: meanN - estimate.toNumber(), se: sdN, interval: { ...toInterval(lo, hi), level: conf }, replicates: reps,
  });
}
export class Bootstrap {
  constructor(fields) {
    Object.assign(this, fields);
    Object.freeze(this);
  }
  toString() {
    const f = (v) => (typeof v === "number" ? String(Number(v.toPrecision(5))) : show(v, 5));
    const how = this.kind === "exact" ? `exact over ${this.distinct} distinct resamples` : `simulated: ${this.resamples} resamples, seed ${JSON.stringify(this.seed)}`;
    return `${this.statistic} ${show(this.estimate, 6)}; bootstrap se ${f(this.se)}, bias ${f(this.bias)}; ${this.method} ${show(this.interval.level, 3)} interval [${show(this.interval.low, 5)}, ${show(this.interval.high, 5)}] (${how})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}
