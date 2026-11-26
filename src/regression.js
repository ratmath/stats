// Regression and analysis of variance, exactly.
//
// Least squares solves the normal equations XᵀX b = Xᵀy in exact rational
// arithmetic, so for rational data the coefficients, fitted values,
// residuals, sums of squares, R², leverages and Cook's distances are exact
// rationals. Standard errors are square roots and so certified; t statistics
// keep their squares exact, F statistics are exact, and every p-value is a
// certified enclosure.
import { refuse, toQ, Q, ZERO, ONE, sqrt, isQ, whole } from "./num.js";
import { TestResult } from "./tests.js";
import { studentT, fDist } from "./continuous.js";
import { exactValues, alternativeOf, sidedP, signedRoot, centredInterval } from "./inference.js";

const sign = (q) => (q.numerator > 0n ? 1 : q.numerator < 0n ? -1 : 0);
const kindOf = (p) => (isQ(p) ? "exact" : "certified");
const dot = (a, b) => a.reduce((s, v, i) => s.add(v.multiply(b[i])), ZERO);

/**
 * The inverse of a symmetric matrix of rationals by Gauss–Jordan
 * elimination. A zero pivot means a column is a combination of earlier ones:
 * its index is reported.
 */
function invert(A) {
  const p = A.length;
  const M = A.map((row, i) => [...row, ...A.map((_, j) => (i === j ? ONE : ZERO))]);
  for (let c = 0; c < p; c++) {
    let r = -1;
    for (let i = c; i < p; i++) if (M[i][c].numerator !== 0n) {
      r = i;
      break;
    }
    if (r < 0) return { singular: c };
    [M[c], M[r]] = [M[r], M[c]];
    const inv = M[c][c].reciprocal();
    M[c] = M[c].map((v) => v.multiply(inv));
    for (let i = 0; i < p; i++) {
      if (i === c || M[i][c].numerator === 0n) continue;
      const f = M[i][c];
      M[i] = M[i].map((v, j) => v.subtract(f.multiply(M[c][j])));
    }
  }
  return { inverse: M.map((row) => row.slice(p)) };
}
/** The first column (in order) that is a combination of the columns before it, or −1. */
function dependentColumn(cols) {
  for (let k = 1; k <= cols.length; k++) {
    const sub = cols.slice(0, k);
    const G = sub.map((a) => sub.map((b) => dot(a, b)));
    if (invert(G).singular !== undefined) return k - 1;
  }
  return -1;
}

/**
 * Columns from what the caller gives:
 *   [1, 2, 3]                 one predictor, named x
 *   [[1, 2], [3, 4], …]       rows of several predictors, named x1, x2, …
 *   { age: […], dose: […] }   named columns
 */
function predictors(x, names) {
  if (x !== null && typeof x === "object" && !Array.isArray(x) && typeof x[Symbol.iterator] !== "function") {
    const keys = Object.keys(x);
    if (!keys.length) refuse("bad-argument", "no predictors were given");
    return { cols: keys.map((k) => exactValues(x[k], `the column ${k}`)), names: keys };
  }
  const list = Array.from(x ?? []);
  if (!list.length) refuse("no-data", "no observations");
  if (Array.isArray(list[0])) {
    const p = list[0].length;
    if (list.some((r) => !Array.isArray(r) || r.length !== p)) refuse("bad-argument", "every row needs the same number of predictors");
    const cols = Array.from({ length: p }, (_, j) => exactValues(list.map((r) => r[j]), `predictor ${j + 1}`));
    return { cols, names: names ?? cols.map((_, j) => `x${j + 1}`) };
  }
  return { cols: [exactValues(list, "x")], names: names ?? ["x"] };
}

/**
 * A fitted linear model y = b₀ + b₁x₁ + … + e.
 *
 *   coefficients, names         the exact least-squares coefficients
 *   fitted, residuals           exact
 *   sse, ssr, sst               residual, regression and total sums of squares
 *   rSquared, adjustedRSquared  exact
 *   df                          residual degrees of freedom n − p
 *   variance                    s² = sse/df, exact; sigma = s, certified
 *   covariance                  s²(XᵀX)⁻¹, exact
 *   leverage, cooksDistance     exact diagnostics, one per observation
 *
 *   table()                      estimate, standard error, t and p for each coefficient
 *   test(name, { value, alternative })    the t test of one coefficient
 *   interval(name, { level })    its confidence interval
 *   fTest()                      the overall F test (all slopes zero)
 *   predict(x)                   the exact fitted value at new predictors
 *   predictionInterval(x, { level, kind })   for the mean ("mean") or a new observation ("new")
 *   standardizedResiduals()      eᵢ/(s√(1 − hᵢ)), certified
 */
export class LinearModel {
  constructor(fields) {
    Object.assign(this, fields);
    Object.freeze(this);
  }
  #index(name) {
    const j = typeof name === "number" ? name : this.names.indexOf(name);
    if (!Number.isInteger(j) || j < 0 || j >= this.names.length) refuse("bad-argument", `no coefficient ${JSON.stringify(name)}; the names are ${this.names.join(", ")}`);
    return j;
  }
  #needsDf() {
    if (this.df < 1) refuse("no-data", "there are no residual degrees of freedom (as many coefficients as observations), so nothing can be inferred");
  }
  /** The exact value of a coefficient by name or position. */
  coefficient(name) {
    return this.coefficients[this.#index(name)];
  }
  /** The t test of H₀: coefficient = value. */
  test(name, { value = 0, alternative = "two-sided", bits = 64 } = {}) {
    this.#needsDf();
    const j = this.#index(name);
    const alt = alternativeOf(alternative);
    const v0 = toQ(value, "the null value");
    const b = this.coefficients[j];
    const se2 = this.covariance[j][j];
    const d = b.subtract(v0);
    if (se2.numerator === 0n) refuse("undefined", "the fit is perfect (no residual variation), so the standard error is 0");
    const t2 = d.multiply(d).divide(se2);
    const pValue = sidedP(fDist(1, this.df).sf(t2, { bits }), sign(d), alt);
    return new TestResult({
      test: `t test of the coefficient ${this.names[j]}`, kind: kindOf(pValue), statistic: signedRoot(t2, sign(d), bits), statisticSquared: t2, df: this.df,
      estimate: b, null: v0, standardError: sqrt(se2, bits), alternative: alt, pValue,
      interval: (options = {}) => this.interval(j, { bits, ...options }),
    });
  }
  /** The confidence interval for one coefficient. */
  interval(name, { level = "0.95", alternative = "two-sided", bits = 64 } = {}) {
    this.#needsDf();
    const j = this.#index(name);
    return centredInterval(this.coefficients[j], this.covariance[j][j], studentT(this.df), { level, alternative, bits, method: "t" });
  }
  /** One row per coefficient: { name, estimate, standardError, t, pValue }. */
  table({ bits = 64 } = {}) {
    return this.names.map((name, j) => {
      const r = this.test(j, { bits });
      return { name, estimate: r.estimate, standardError: r.standardError, t: r.statistic, pValue: r.pValue };
    });
  }
  /** The overall F test: do the predictors explain anything beyond the intercept? */
  fTest({ bits = 64 } = {}) {
    this.#needsDf();
    const k = this.intercept ? this.p - 1 : this.p;
    if (k < 1) refuse("bad-argument", "a model with only an intercept has no slopes to test");
    if (this.sse.numerator === 0n) return new TestResult({ test: "F test of the regression", kind: "exact", statistic: null, df: [k, this.df], alternative: "greater", pValue: ZERO, reason: "the fit is perfect" });
    const F = this.ssr.divide(Q(k)).divide(this.variance);
    const pValue = fDist(k, this.df).sf(F, { bits });
    return new TestResult({ test: "F test of the regression", kind: kindOf(pValue), statistic: F, df: [k, this.df], alternative: "greater", pValue });
  }
  #row(x) {
    let raw;
    const k = this.predictorNames.length;
    if (x !== null && typeof x === "object" && !Array.isArray(x)) raw = this.predictorNames.map((n) => {
      if (x[n] === undefined) refuse("bad-argument", `the new point has no value for ${n}`);
      return x[n];
    });
    else raw = Array.isArray(x) ? x : [x];
    if (raw.length !== k) refuse("bad-argument", `the model has ${k} predictor${k === 1 ? "" : "s"}, the point has ${raw.length}`);
    const vals = raw.map((v) => toQ(v, "a predictor value"));
    const row = this.expand(vals);
    return this.intercept ? [ONE, ...row] : row;
  }
  /** The fitted value at new predictor values, exactly. */
  predict(x) {
    return dot(this.coefficients, this.#row(x));
  }
  /**
   * An interval at new predictor values: for the mean response (kind
   * "mean", the default) or for a single new observation (kind "new", which
   * adds s² for its own error).
   */
  predictionInterval(x, { level = "0.95", kind = "mean", bits = 64 } = {}) {
    this.#needsDf();
    if (kind !== "mean" && kind !== "new") refuse("bad-argument", `kind is "mean" or "new", got "${kind}"`);
    const row = this.#row(x);
    const h = dot(row, this.inverse.map((r) => dot(r, row)));
    const se2 = this.variance.multiply(kind === "new" ? h.add(ONE) : h);
    return centredInterval(dot(this.coefficients, row), se2, studentT(this.df), { level, bits, method: kind === "new" ? "prediction" : "confidence" });
  }
  /** eᵢ/(s√(1 − hᵢ)): certified, null where hᵢ = 1. */
  standardizedResiduals({ bits = 64 } = {}) {
    this.#needsDf();
    return this.residuals.map((e, i) => {
      const one = ONE.subtract(this.leverage[i]);
      if (one.numerator === 0n || this.variance.numerator === 0n) return null;
      const s = e.multiply(e).divide(this.variance.multiply(one));
      return signedRoot(s, sign(e), bits);
    });
  }
  get sigma() {
    return this.df < 1 ? null : sqrt(this.variance, 64);
  }
  toString() {
    const terms = this.coefficients.map((b, j) => (this.names[j] === "(intercept)" ? `${b}` : `${b}·${this.names[j]}`));
    return `y = ${terms.join(" + ").replace(/\+ -/g, "− ")}  (n = ${this.n}, R² = ${this.rSquared})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/**
 * Fit y = b₀ + Σ bⱼxⱼ by least squares, exactly.
 *   linearModel(x, y)                         x: a list, rows, or { name: column }
 *   linearModel(x, y, { intercept: false })   through the origin
 * Refuses, naming the column, when a predictor is a linear combination of
 * the others (the coefficients would not be determined).
 */
export function linearModel(x, y, { intercept = true, names, expand, predictorNames } = {}) {
  const ys = exactValues(y, "y");
  const pr = predictors(x, names);
  const n = ys.length;
  if (pr.cols.some((c) => c.length !== n)) refuse("bad-argument", `y has ${n} values but a predictor has ${pr.cols.find((c) => c.length !== n).length}`);
  const cols = intercept ? [ys.map(() => ONE), ...pr.cols] : pr.cols;
  const allNames = intercept ? ["(intercept)", ...pr.names] : pr.names;
  const p = cols.length;
  if (n < p) refuse("no-data", `${p} coefficients need at least ${p} observations, there are ${n}`);
  const XtX = cols.map((a) => cols.map((b) => dot(a, b)));
  const inv = invert(XtX);
  if (inv.inverse === undefined) {
    const j = dependentColumn(cols);
    refuse("collinear", `${allNames[j]} is a linear combination of ${j === 0 ? "nothing (it is all zero)" : allNames.slice(0, j).join(", ")}, so the coefficients are not determined`);
  }
  const A = inv.inverse;
  const Xty = cols.map((c) => dot(c, ys));
  const coefficients = A.map((row) => dot(row, Xty));
  const rows = ys.map((_, i) => cols.map((c) => c[i]));
  const fitted = rows.map((r) => dot(r, coefficients));
  const residuals = ys.map((v, i) => v.subtract(fitted[i]));
  const sse = residuals.reduce((s, e) => s.add(e.multiply(e)), ZERO);
  const ybar = ys.reduce((s, v) => s.add(v), ZERO).divide(Q(n));
  const sst = intercept ? ys.reduce((s, v) => s.add(v.subtract(ybar).multiply(v.subtract(ybar))), ZERO) : ys.reduce((s, v) => s.add(v.multiply(v)), ZERO);
  const ssr = sst.subtract(sse);
  const df = n - p;
  const variance = df > 0 ? sse.divide(Q(df)) : null;
  const leverage = rows.map((r) => dot(r, A.map((row) => dot(row, r))));
  const rSquared = sst.numerator === 0n ? ONE : ONE.subtract(sse.divide(sst));
  const denomDf = intercept ? n - 1 : n;
  const adjustedRSquared = df > 0 && sst.numerator !== 0n ? ONE.subtract(sse.divide(Q(df)).divide(sst.divide(Q(denomDf)))) : null;
  const cooksDistance = rows.map((_, i) => {
    if (!variance || variance.numerator === 0n) return null;
    const one = ONE.subtract(leverage[i]);
    if (one.numerator === 0n) return null;
    const e = residuals[i];
    return e.multiply(e).multiply(leverage[i]).divide(Q(p).multiply(variance).multiply(one).multiply(one));
  });
  return new LinearModel({
    n, p, df, intercept, names: allNames, predictorNames: predictorNames ?? pr.names, expand: expand ?? ((v) => v),
    coefficients, fitted, residuals, sse, ssr, sst, rSquared, adjustedRSquared, variance,
    covariance: variance ? A.map((row) => row.map((v) => v.multiply(variance))) : null,
    inverse: A, leverage, cooksDistance, y: ys, columns: cols,
  });
}

/**
 * Least-squares polynomial of the given degree: y = b₀ + b₁x + … + b_d·x^d,
 * exactly. predict takes a single x.
 */
export function polyFit(x, y, degree, options = {}) {
  const d = whole(degree, "the degree", { min: 1 });
  const xs = exactValues(x, "x");
  const powers = (v) => Array.from({ length: d }, (_, k) => v.pow(k + 1));
  const names = Array.from({ length: d }, (_, k) => (k === 0 ? "x" : `x^${k + 1}`));
  return linearModel(xs.map(powers), y, { ...options, names, expand: (vals) => powers(vals[0]), predictorNames: ["x"] });
}

/**
 * 0/1 indicator columns for a categorical variable, one for each level but
 * the reference (the first level in order of appearance, unless given), for
 * use as predictors: linearModel({ ...dummies(group), dose }, y).
 */
export function dummies(labels, { reference, name = "" } = {}) {
  const ls = Array.from(labels, String);
  if (!ls.length) refuse("no-data", "no labels");
  const levels = [...new Set(ls)];
  const ref = reference === undefined ? levels[0] : String(reference);
  if (!levels.includes(ref)) refuse("bad-argument", `the reference level ${ref} does not occur`);
  const out = {};
  for (const l of levels) if (l !== ref) out[`${name}${name ? "=" : ""}${l}`] = ls.map((v) => (v === l ? 1 : 0));
  return out;
}

/**
 * The F test of a smaller model against a larger one that contains it (the
 * same y, and every column of the smaller model in the column space of the
 * larger, which is checked exactly):
 * F = ((sse₀ − sse₁)/(p₁ − p₀)) / (sse₁/df₁).
 */
export function compareModels(small, big, { bits = 64 } = {}) {
  if (!(small instanceof LinearModel) || !(big instanceof LinearModel)) refuse("bad-argument", "compareModels needs two linear models");
  if (small.n !== big.n || small.y.some((v, i) => !v.equals(big.y[i]))) refuse("bad-argument", "the models must be fitted to the same y");
  if (small.p >= big.p) refuse("bad-argument", "the first model must have fewer coefficients than the second");
  // every column of the smaller model must be reproduced exactly by the larger one's columns
  for (const [j, col] of small.columns.entries()) {
    const fit = big.inverse.map((row) => dot(row, big.columns.map((c) => dot(c, col))));
    if (col.some((v, i) => !v.equals(dot(big.columns.map((c) => c[i]), fit)))) refuse("not-nested", `the column ${small.names[j]} of the smaller model is not in the span of the larger model's columns`);
  }
  if (big.df < 1) refuse("no-data", "the larger model has no residual degrees of freedom");
  const k = big.p - small.p;
  const drop = small.sse.subtract(big.sse);
  if (big.sse.numerator === 0n) return new TestResult({ test: "F test of nested models", kind: "exact", statistic: null, df: [k, big.df], alternative: "greater", pValue: drop.numerator === 0n ? ONE : ZERO, reason: "the larger model fits perfectly" });
  const F = drop.divide(Q(k)).divide(big.variance);
  const pValue = fDist(k, big.df).sf(F, { bits });
  return new TestResult({ test: "F test of nested models", kind: kindOf(pValue), statistic: F, df: [k, big.df], alternative: "greater", pValue, sseDrop: drop });
}

/**
 * One-way analysis of variance: do the groups share a mean?
 *   anova([[…], […], …])  or  anova({ a: […], b: […] })
 * The sums of squares, mean squares and F are exact; the result also has
 * `table` (rows for between, within and total), the group means, and
 * etaSquared = SSB/SST.
 */
export function anova(groups, { bits = 64 } = {}) {
  const named = groups !== null && typeof groups === "object" && !Array.isArray(groups);
  const labels = named ? Object.keys(groups) : Array.from(groups ?? [], (_, i) => `group ${i + 1}`);
  const gs = (named ? Object.values(groups) : Array.from(groups ?? [])).map((g, i) => exactValues(g, labels[i]));
  const k = gs.length;
  if (k < 2) refuse("bad-argument", "analysis of variance needs at least two groups");
  const N = gs.reduce((s, g) => s + g.length, 0);
  if (N <= k) refuse("no-data", "there must be more observations than groups");
  const means = gs.map((g) => g.reduce((s, v) => s.add(v), ZERO).divide(Q(g.length)));
  const grand = gs.flat().reduce((s, v) => s.add(v), ZERO).divide(Q(N));
  const ssb = gs.reduce((s, g, i) => s.add(Q(g.length).multiply(means[i].subtract(grand).pow(2))), ZERO);
  const ssw = gs.reduce((s, g, i) => s.add(g.reduce((t, v) => t.add(v.subtract(means[i]).pow(2)), ZERO)), ZERO);
  const sst = ssb.add(ssw);
  const [df1, df2] = [k - 1, N - k];
  const table = [
    { source: "between groups", df: df1, ss: ssb, ms: ssb.divide(Q(df1)) },
    { source: "within groups", df: df2, ss: ssw, ms: ssw.divide(Q(df2)) },
    { source: "total", df: N - 1, ss: sst, ms: null },
  ];
  const extra = { table, means: Object.fromEntries(labels.map((l, i) => [l, means[i]])), grandMean: grand, etaSquared: sst.numerator === 0n ? null : ssb.divide(sst) };
  if (ssw.numerator === 0n) {
    return new TestResult({ test: "one-way analysis of variance", kind: "exact", statistic: null, df: [df1, df2], alternative: "greater", pValue: ssb.numerator === 0n ? ONE : ZERO, reason: "there is no variation within the groups", ...extra });
  }
  const F = table[0].ms.divide(table[1].ms);
  const pValue = fDist(df1, df2).sf(F, { bits });
  return new TestResult({ test: "one-way analysis of variance", kind: kindOf(pValue), statistic: F, df: [df1, df2], alternative: "greater", pValue, ...extra });
}
