// Descriptive statistics, exactly. Data are rationals (JS numbers are read by
// their decimal form: 0.1 is 1/10) or intervals for measurements with stated
// uncertainty; with intervals every answer is an enclosure of the statistic
// over all the values the data could have.
import { RationalInterval } from "@ratmath/core";
import { refuse, toQ, toValue, ZERO, ONE, Q, add, sub, mul, div, square, sqrt, lowOf, highOf, interval, isQ, whole, probability, min as minQ, max as maxQ } from "./num.js";

/** The data as exact values, refusing an empty list. */
export function values(data, what = "the data") {
  if (data == null || typeof data[Symbol.iterator] !== "function" || typeof data === "string") refuse("bad-argument", `${what} must be a list of numbers`);
  const out = Array.from(data, (v) => toValue(v, "a data value"));
  if (!out.length) refuse("no-data", `${what} is empty`);
  return out;
}
const exactOnly = (xs, what) => {
  if (xs.some((x) => !isQ(x))) refuse("needs-exact", `${what} is not defined for interval data`);
  return xs;
};
const sorted = (xs) => [...xs].sort((a, b) => a.compareTo(b));
/** A statistic that increases with every data value, on interval data: its value at the low ends and at the high ends. */
const monotone = (xs, f) => {
  if (xs.every(isQ)) return f(sorted(xs));
  const lo = f(sorted(xs.map(lowOf)));
  const hi = f(sorted(xs.map(highOf)));
  return lo.equals(hi) ? lo : interval(lo, hi);
};

export const count = (data) => values(data).length;
export const sum = (data) => values(data).reduce(add, ZERO);
export const mean = (data) => {
  const xs = values(data);
  return div(xs.reduce(add, ZERO), Q(xs.length));
};
/** Σ wᵢxᵢ / Σ wᵢ. */
export function weightedMean(data, weights) {
  const xs = values(data);
  const ws = values(weights, "the weights");
  if (ws.length !== xs.length) refuse("bad-argument", "one weight for each value");
  return div(xs.reduce((s, x, i) => add(s, mul(x, ws[i])), ZERO), ws.reduce(add, ZERO));
}
/**
 * The variance: the sum of squared deviations from the mean over n − 1, or
 * over n with { population: true }.
 */
export function variance(data, { population = false } = {}) {
  const xs = values(data);
  const n = xs.length;
  if (!population && n < 2) refuse("no-data", "the sample variance needs at least two values");
  const m = div(xs.reduce(add, ZERO), Q(n));
  const ss = xs.reduce((s, x) => add(s, square(sub(x, m))), ZERO);
  return div(ss, Q(population ? n : n - 1));
}
/** The standard deviation: exact when the variance is a perfect square, otherwise certified. */
export const sd = (data, { population = false, bits = 96 } = {}) => sqrt(variance(data, { population }), bits);
/** The k-th moment about 0, or about the mean with { central: true } (divisor n). */
export function moment(data, k, { central = false } = {}) {
  const xs = exactOnly(values(data), "moment");
  const p = whole(k, "the order of the moment");
  const m = central ? mean(xs) : ZERO;
  return xs.reduce((s, x) => s.add(x.subtract(m).pow(p)), ZERO).divide(Q(xs.length));
}
/** The skewness m₃/m₂^(3/2) (certified unless it is rational). */
export function skewness(data, { bits = 96 } = {}) {
  const m2 = moment(data, 2, { central: true });
  if (m2.numerator === 0n) refuse("undefined", "the skewness is not defined when all the values are equal");
  const m3 = moment(data, 3, { central: true });
  const r = sqrt(m3.multiply(m3).divide(m2.pow(3)), bits);
  if (m3.numerator >= 0n) return r;
  return isQ(r) ? r.negate() : new RationalInterval(r.high.negate(), r.low.negate());
}
/** The kurtosis m₄/m₂², exactly; 3 for a normal distribution. With { excess: true }, minus 3. */
export function kurtosis(data, { excess = false } = {}) {
  const m2 = moment(data, 2, { central: true });
  if (m2.numerator === 0n) refuse("undefined", "the kurtosis is not defined when all the values are equal");
  const k = moment(data, 4, { central: true }).divide(m2.multiply(m2));
  return excess ? k.subtract(Q(3)) : k;
}

export const min = (data) => monotone(values(data), (s) => s[0]);
export const max = (data) => monotone(values(data), (s) => s[s.length - 1]);
export const range = (data) => {
  const xs = values(data);
  if (xs.every(isQ)) return sub(max(xs), min(xs));
  // the largest and smallest possible spread
  const hi = sub(highOf(max(xs)), lowOf(min(xs)));
  const top = xs.map(lowOf).reduce(maxQ);
  const bottom = xs.map(highOf).reduce(minQ);
  const lo = maxQ(ZERO, top.subtract(bottom));
  return interval(lo, hi);
};

// name: [a, b] in the rule h = (n + a)·p + b, counting positions from 1 (Hyndman and Fan's types)
const METHODS = {
  linear: [Q(-1), ONE], // type 7, the default of most software
  weibull: [ONE, ZERO], // type 6
  hazen: [ZERO, Q(1, 2)], // type 5
  "median-unbiased": [Q(1, 3), Q(1, 3)], // type 8
};
function quantileSorted(s, p, method) {
  const n = s.length;
  if (method === "lower" || method === "higher" || method === "nearest" || method === "inverse") {
    // by counting: the smallest value with at least a fraction p of the data at or below it
    const np = p.multiply(Q(n));
    const fl = Number(np.floor());
    const whole0 = np.denominator === 1n;
    if (method === "inverse") return s[Math.min(n - 1, Math.max(0, whole0 ? fl - 1 : fl))];
    const h = p.multiply(Q(n - 1));
    const lo = Number(h.floor());
    const frac = h.subtract(Q(lo));
    if (method === "lower" || frac.numerator === 0n) return s[lo];
    if (method === "higher") return s[lo + 1];
    return frac.compareTo(Q(1, 2)) <= 0 ? s[lo] : s[lo + 1];
  }
  if (method === "midpoint") {
    const h = p.multiply(Q(n - 1));
    const lo = Number(h.floor());
    return h.denominator === 1n ? s[lo] : s[lo].add(s[lo + 1]).divide(Q(2));
  }
  const rule = METHODS[method];
  if (!rule) refuse("bad-argument", `unknown quantile method "${method}"; use linear, weibull, hazen, median-unbiased, inverse, lower, higher, nearest or midpoint`);
  const [a, b] = rule;
  const h = Q(n).add(a).multiply(p).add(b); // a position from 1 to n
  if (h.compareTo(ONE) <= 0) return s[0];
  if (h.compareTo(Q(n)) >= 0) return s[n - 1];
  const lo = Number(h.floor());
  const frac = h.subtract(Q(lo));
  return s[lo - 1].add(s[lo].subtract(s[lo - 1]).multiply(frac));
}
/**
 * The p-quantile of the data. The default method "linear" interpolates
 * between the two nearest values in sorted order (type 7); others are
 * "weibull", "hazen", "median-unbiased", "inverse" (the empirical
 * distribution's own quantile), "lower", "higher", "nearest" and "midpoint".
 */
export function quantile(data, p, { method = "linear" } = {}) {
  const level = probability(p, "the level");
  return monotone(values(data), (s) => quantileSorted(s, level, method));
}
export const median = (data) => quantile(data, Q(1, 2));
/** [Q1, median, Q3]. */
export const quartiles = (data, options) => [Q(1, 4), Q(1, 2), Q(3, 4)].map((p) => quantile(data, p, options));
/** Q3 − Q1 (for interval data, an enclosure). */
export const iqr = (data, options) => {
  const [q1, , q3] = quartiles(data, options);
  return sub(q3, q1);
};
/** [minimum, Q1, median, Q3, maximum]. */
export const fiveNumber = (data, options) => {
  const [q1, q2, q3] = quartiles(data, options);
  return [min(data), q1, q2, q3, max(data)];
};
/** The median of the absolute deviations from the median. */
export function mad(data) {
  const xs = exactOnly(values(data), "mad");
  const m = median(xs);
  return median(xs.map((x) => x.subtract(m).abs()));
}
/** The mean after dropping the given fraction of values from each end (rounded down to whole values). */
export function trimmedMean(data, fraction) {
  const xs = values(data);
  const f = toQ(fraction, "the fraction trimmed");
  if (f.numerator < 0n || f.compareTo(Q(1, 2)) >= 0) refuse("bad-argument", "the fraction trimmed from each end must be at least 0 and below 1/2");
  const k = Number(f.multiply(Q(xs.length)).floor());
  return monotone(xs, (s) => s.slice(k, s.length - k).reduce((a, b) => a.add(b), ZERO).divide(Q(s.length - 2 * k)));
}

/** Rows { value, count, proportion, cumulative } for each distinct value, in increasing order. */
export function frequencies(data) {
  const xs = sorted(exactOnly(values(data), "frequencies"));
  const n = Q(xs.length);
  const rows = [];
  for (const x of xs) {
    const last = rows[rows.length - 1];
    if (last && last.value.equals(x)) last.count++;
    else rows.push({ value: x, count: 1 });
  }
  let running = 0;
  for (const r of rows) {
    running += r.count;
    r.proportion = Q(r.count).divide(n);
    r.cumulative = Q(running).divide(n);
  }
  return rows;
}
/** The most frequent values (all of them when there is a tie). */
export function modes(data) {
  const rows = frequencies(data);
  const best = Math.max(...rows.map((r) => r.count));
  return rows.filter((r) => r.count === best).map((r) => r.value);
}
/**
 * Counts in classes [edge, next edge), the last class closed on the right.
 * Give { edges }, or { width } (and optionally { start }), or { bins }.
 */
export function histogram(data, { edges, width, start, bins } = {}) {
  const xs = sorted(exactOnly(values(data), "histogram"));
  const lo = xs[0];
  const hi = xs[xs.length - 1];
  let cut;
  if (edges) cut = Array.from(edges, (e) => toQ(e, "an edge"));
  else {
    const first = start !== undefined ? toQ(start, "start") : lo;
    let w;
    if (width !== undefined) w = toQ(width, "width");
    else {
      const k = bins !== undefined ? whole(bins, "bins", { min: 1 }) : Math.max(1, Math.ceil(Math.log2(xs.length) + 1));
      w = hi.subtract(first).divide(Q(k));
      if (w.numerator === 0n) w = ONE;
    }
    if (w.numerator <= 0n) refuse("bad-argument", "the class width must be positive");
    cut = [first];
    while (cut[cut.length - 1].compareTo(hi) < 0) cut.push(cut[cut.length - 1].add(w));
    if (cut.length === 1) cut.push(first.add(w));
  }
  if (cut.length < 2) refuse("bad-argument", "a histogram needs at least two edges");
  const rows = [];
  for (let i = 0; i + 1 < cut.length; i++) {
    if (cut[i + 1].compareTo(cut[i]) <= 0) refuse("bad-argument", "the edges must increase");
    rows.push({ from: cut[i], to: cut[i + 1], count: 0 });
  }
  let outside = 0;
  for (const x of xs) {
    let placed = false;
    for (let i = 0; i < rows.length && !placed; i++) {
      const last = i === rows.length - 1;
      if (x.compareTo(rows[i].from) >= 0 && (x.compareTo(rows[i].to) < 0 || (last && x.equals(rows[i].to)))) {
        rows[i].count++;
        placed = true;
      }
    }
    if (!placed) outside++;
  }
  const n = Q(xs.length);
  for (const r of rows) {
    r.proportion = Q(r.count).divide(n);
    r.density = r.proportion.divide(r.to.subtract(r.from));
  }
  return Object.assign(rows, { outside });
}

const pairs = (x, y) => {
  const xs = values(x, "x");
  const ys = values(y, "y");
  if (xs.length !== ys.length) refuse("bad-argument", `x has ${xs.length} values and y has ${ys.length}`);
  return [xs, ys];
};
/** The covariance, with divisor n − 1 (or n with { population: true }). */
export function covariance(x, y, { population = false } = {}) {
  const [xs, ys] = pairs(x, y);
  const n = xs.length;
  if (!population && n < 2) refuse("no-data", "the sample covariance needs at least two pairs");
  const mx = mean(xs);
  const my = mean(ys);
  return div(xs.reduce((s, v, i) => add(s, mul(sub(v, mx), sub(ys[i], my))), ZERO), Q(population ? n : n - 1));
}
/** Pearson's correlation coefficient: exact when it is rational, otherwise certified. */
export function correlation(x, y, { bits = 96 } = {}) {
  const [xs, ys] = pairs(x, y);
  exactOnly(xs, "correlation");
  exactOnly(ys, "correlation");
  const c = covariance(xs, ys, { population: true });
  const vv = variance(xs, { population: true }).multiply(variance(ys, { population: true }));
  if (vv.numerator === 0n) refuse("undefined", "the correlation is not defined when one variable is constant");
  const r = sqrt(c.multiply(c).divide(vv), bits);
  if (c.numerator >= 0n) return r;
  return isQ(r) ? r.negate() : new RationalInterval(r.high.negate(), r.low.negate());
}
/**
 * The least-squares line y = intercept + slope·x, exactly.
 *   { slope, intercept, r2, residuals, fitted, sse, sst, n, predict(x) }
 */
export function leastSquares(x, y) {
  const [xs, ys] = pairs(x, y);
  exactOnly(xs, "leastSquares");
  exactOnly(ys, "leastSquares");
  const n = xs.length;
  if (n < 2) refuse("no-data", "a line needs at least two points");
  const mx = mean(xs);
  const my = mean(ys);
  let sxx = ZERO;
  let sxy = ZERO;
  let syy = ZERO;
  for (let i = 0; i < n; i++) {
    const dx = xs[i].subtract(mx);
    const dy = ys[i].subtract(my);
    sxx = sxx.add(dx.multiply(dx));
    sxy = sxy.add(dx.multiply(dy));
    syy = syy.add(dy.multiply(dy));
  }
  if (sxx.numerator === 0n) refuse("undefined", "all the x values are equal, so no line is determined");
  const slope = sxy.divide(sxx);
  const intercept = my.subtract(slope.multiply(mx));
  const predict = (v) => intercept.add(slope.multiply(toQ(v, "x")));
  const fitted = xs.map(predict);
  const residuals = ys.map((v, i) => v.subtract(fitted[i]));
  const sse = residuals.reduce((s, r) => s.add(r.multiply(r)), ZERO);
  return { slope, intercept, r2: syy.numerator === 0n ? ONE : ONE.subtract(sse.divide(syy)), residuals, fitted, sse, sst: syy, n, predict };
}

/**
 * Everything at once: n, mean, variance, sd, the five-number summary, the
 * interquartile range and the values beyond 1.5 of it from the quartiles.
 */
export function summary(data, { bits = 96 } = {}) {
  const xs = values(data);
  const [lo, q1, med, q3, hi] = fiveNumber(xs);
  const out = { n: xs.length, mean: mean(xs), min: lo, q1, median: med, q3, max: hi, iqr: sub(q3, q1), range: range(xs) };
  if (xs.length > 1) {
    out.variance = variance(xs);
    out.sd = sqrt(out.variance, bits);
  }
  if (xs.every(isQ)) {
    const fence = out.iqr.multiply(Q(3, 2));
    out.outliers = sorted(xs).filter((x) => x.compareTo(q1.subtract(fence)) < 0 || x.compareTo(q3.add(fence)) > 0);
    out.modes = modes(xs);
  }
  return out;
}

// ---- the same statistics on plain JS numbers, for simulations where speed matters ----

const fsum = (a) => {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i];
  return s;
};
const fsorted = (a) => Float64Array.from(a).sort();
const fquantile = (a, p) => {
  const s = fsorted(a);
  const h = (s.length - 1) * p;
  const lo = Math.floor(h);
  return lo + 1 < s.length ? s[lo] + (h - lo) * (s[lo + 1] - s[lo]) : s[lo];
};
/** Floating-point versions for simulation: quick, approximate, and never mixed with the exact ones. */
export const fast = Object.freeze({
  sum: fsum,
  mean: (a) => fsum(a) / a.length,
  variance(a, { population = false } = {}) {
    const m = fsum(a) / a.length;
    let s = 0;
    for (let i = 0; i < a.length; i++) s += (a[i] - m) ** 2;
    return s / (population ? a.length : a.length - 1);
  },
  sd: (a, options) => Math.sqrt(fast.variance(a, options)),
  quantile: fquantile,
  median: (a) => fquantile(a, 0.5),
  min: (a) => Math.min(...a),
  max: (a) => Math.max(...a),
  trimmedMean(a, fraction) {
    const s = fsorted(a);
    const k = Math.floor(fraction * s.length);
    return fsum(s.subarray(k, s.length - k)) / (s.length - 2 * k);
  },
  covariance(x, y) {
    const mx = fsum(x) / x.length;
    const my = fsum(y) / y.length;
    let s = 0;
    for (let i = 0; i < x.length; i++) s += (x[i] - mx) * (y[i] - my);
    return s / (x.length - 1);
  },
  correlation: (x, y) => fast.covariance(x, y) / Math.sqrt(fast.variance(x) * fast.variance(y)),
  leastSquares(x, y) {
    const mx = fsum(x) / x.length;
    const my = fsum(y) / y.length;
    let sxx = 0;
    let sxy = 0;
    let syy = 0;
    for (let i = 0; i < x.length; i++) {
      sxx += (x[i] - mx) ** 2;
      sxy += (x[i] - mx) * (y[i] - my);
      syy += (y[i] - my) ** 2;
    }
    const slope = sxy / sxx;
    const intercept = my - slope * mx;
    return { slope, intercept, r2: syy === 0 ? 1 : (sxy * sxy) / (sxx * syy), n: x.length, predict: (v) => intercept + slope * v };
  },
});
