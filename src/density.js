// Continuous distributions from a density the user writes, and new
// distributions from old ones: a monotone function of a random variable and
// the order statistics of a sample.
//
// A density is first proven to be one: f ≥ 0 on the support and total mass 1.
// For a piecewise polynomial both are decided exactly (the sign between the
// isolated roots, the mass from the antiderivative), and the CDF, the moments
// and the order statistics stay exact. For any other formula the values come
// from analytic's certified integrals, and the mass must be enclosed around 1.
import { Rational, RationalInterval } from "@ratmath/core";
import { parse, format, fn, zeros, encloseAst, differentiate, simplify, freeVariables, compileFloat, exactOf } from "@ratmath/analytic";
import { Continuous, u01 } from "./continuous.js";
import { crossing } from "./solve.js";
import { refuse, toQ, ZERO, ONE, Q, add, sub, mul, div, lowOf, highOf, interval, isQ } from "./num.js";
import * as F from "./float.js";

let IvClass = null;
const ivOf = (lo, hi = lo) => {
  IvClass ??= encloseAst(parse("0"), {}).constructor;
  return new IvClass(lo, hi);
};
const tidyIv = (v) => (v instanceof RationalInterval && v.low.equals(v.high) ? v.low : v);
const astOf = (f) => (typeof f === "string" ? parse(f) : typeof f === "number" ? parse(String(f)) : f);
const endOf = (v, what) => {
  if (v === null || v === "inf" || v === "-inf" || v === Infinity || v === -Infinity) return null;
  return toQ(v, what);
};

// ---- exact polynomials (ascending Rational coefficients) -------------------------------------

const pEval = (p, x) => p.reduceRight((acc, c) => acc.multiply(x).add(c), ZERO);
const pAdd = (p, q) => Array.from({ length: Math.max(p.length, q.length) }, (_, i) => (p[i] ?? ZERO).add(q[i] ?? ZERO));
const pScale = (p, c) => p.map((a) => a.multiply(c));
const pMul = (p, q) => {
  const out = Array(p.length + q.length - 1).fill(ZERO);
  p.forEach((a, i) => q.forEach((b, j) => (out[i + j] = out[i + j].add(a.multiply(b)))));
  return out;
};
const pPow = (p, k) => {
  let out = [ONE];
  for (let i = 0; i < k; i++) out = pMul(out, p);
  return out;
};
/** the antiderivative that vanishes at 0 */
const pInt = (p) => [ZERO, ...p.map((c, i) => c.divide(Q(i + 1)))];
const pIntegral = (p, a, b) => {
  const P = pInt(p);
  return pEval(P, b).subtract(pEval(P, a));
};
const pFloat = (p) => {
  const cs = p.map((c) => c.toNumber());
  return (x) => cs.reduceRight((acc, c) => acc * x + c, 0);
};
const pText = (p) => {
  const parts = [];
  p.forEach((c, k) => {
    if (c.numerator === 0n) return;
    const m = k === 0 ? "" : k === 1 ? "x" : `x^${k}`;
    parts.push(m === "" ? `${c}` : c.equals(ONE) ? m : `${c}*${m}`);
  });
  return parts.length ? parts.join(" + ").replaceAll("+ -", "- ") : "0";
};

/** The coefficients of a polynomial tree in x, exactly, or null when it is not a polynomial. */
function polynomialOf(ast, x) {
  const vars = [...freeVariables(ast)];
  if (vars.some((v) => v !== x)) refuse("bad-argument", `the density may only use the variable ${x}, not ${vars.filter((v) => v !== x).join(", ")}`);
  const out = [];
  let d = ast;
  let fact = ONE;
  for (let k = 0; k <= 64; k++) {
    let c;
    try {
      c = encloseAst(d, { [x]: ivOf(ZERO) }, 64);
    } catch {
      return null;
    }
    if (!c.isPoint) return null;
    out.push(c.lo.divide(fact));
    d = simplify(differentiate(d, x));
    fact = fact.multiply(Q(k + 1));
    if (d.type === "num" && d.value.numerator === 0n) return out;
  }
  return null;
}

/** f ≥ 0 on [a, b] for a polynomial, decided exactly: the sign at points between its isolated roots. */
function polynomialNonnegative(p, ast, a, b, x) {
  const roots = zeros(ast, [a, b], { variable: x }).exact().filter((z) => z.lo !== undefined);
  const samples = [a, b];
  const ends = [a];
  for (const z of roots) ends.push(toQ(z.lo), toQ(z.hi));
  ends.push(b);
  for (let i = 0; i + 1 < ends.length; i += 2) if (ends[i].compareTo(ends[i + 1]) < 0) samples.push(ends[i].add(ends[i + 1]).divide(Q(2)));
  for (const s of samples) if (pEval(p, s).numerator < 0n) return { ok: false, at: s, value: pEval(p, s) };
  return { ok: true };
}

// ---- building a distribution from pieces ---------------------------------------------------

/**
 * A distribution from a density.
 *
 *   fromDensity("3/8*x^2", [0, 2])
 *   fromDensity("x^2", [0, 2], { normalize: true })        // finds c = 3/8
 *   fromDensity("exp(-x)", [0, "inf"])
 *   fromDensity([["x", [0, 1]], ["2 - x", [1, 2]]])         // a piecewise density
 *
 * The density is proven first: f ≥ 0 on the support and total mass 1
 * (exactly for piecewise polynomials, within a certified enclosure
 * otherwise), and refused if not. With { normalize: true } the formula is
 * divided by its mass, which must then be exact; `constant` is the factor.
 * The result is a Continuous distribution: pdf, cdf, quantile, mean,
 * variance, samples, and `transform` and `orderStatistic` like any other.
 * `checks` says how f ≥ 0 and the mass were established.
 */
export function fromDensity(f, support, options = {}) {
  const x = options.variable ?? "x";
  const list = Array.isArray(f) && Array.isArray(f[0]) ? f : [[f, support]];
  if (!list.length) refuse("bad-argument", "fromDensity needs a density");
  let pieces = list.map(([g, [a, b]], i) => {
    const ast = astOf(g);
    const lo = endOf(a, "the start of the support");
    const hi = endOf(b, "the end of the support");
    if (lo !== null && hi !== null && hi.compareTo(lo) <= 0) refuse("bad-argument", `piece ${i + 1} has an empty support`);
    return { ast, a: lo, b: hi, poly: polynomialOf(ast, x) };
  });
  for (let i = 1; i < pieces.length; i++) {
    const [p, q] = [pieces[i - 1], pieces[i]];
    if (p.b === null || q.a === null || !p.b.equals(q.a)) refuse("bad-argument", "the pieces must be listed in order and meet end to end");
  }
  // f ≥ 0
  const checks = { nonnegative: "proven", mass: "exact" };
  for (const p of pieces) {
    if (p.poly) {
      if (p.a === null || p.b === null) refuse("not-a-density", `the polynomial ${format(p.ast)} has infinite mass on an unbounded support`);
      const r = polynomialNonnegative(p.poly, p.ast, p.a, p.b, x);
      if (!r.ok) refuse("not-a-density", `the density is negative at ${x} = ${r.at}: ${r.value}`);
      continue;
    }
    const lo = p.a ?? Q(-64);
    const hi = p.b ?? (p.a ?? ZERO).add(Q(64));
    const r = nonnegativeOn(p.ast, x, lo, hi);
    if (r === "negative") refuse("not-a-density", `the density ${format(p.ast)} is negative somewhere on [${lo}, ${hi}]`);
    if (r !== true) refuse("not-a-density", `could not prove ${format(p.ast)} ≥ 0 on [${lo}, ${hi}]`);
    if (p.a === null || p.b === null) checks.nonnegative = `proven on [${lo}, ${hi}]; beyond that by the sign of the leading term only`;
  }
  // the mass
  let mass = ZERO;
  for (const p of pieces) {
    if (p.poly) mass = add(mass, pIntegral(p.poly, p.a, p.b));
    else {
      mass = add(mass, integralOf(p.ast, x, p.a, p.b));
      checks.mass = "enclosure";
    }
  }
  let constant = ONE;
  if (options.normalize) {
    if (!isQ(mass)) refuse("not-exact", "normalize needs an exact mass, which only a piecewise polynomial has");
    if (mass.numerator <= 0n) refuse("not-a-density", "the total mass is 0");
    constant = mass.reciprocal();
    pieces = pieces.map((p) => ({ ...p, ast: { type: "op", op: "*", args: [{ type: "num", value: constant }, p.ast] }, poly: p.poly && pScale(p.poly, constant) }));
  } else if (isQ(mass) ? !mass.equals(ONE) : !(mass.low.compareTo(ONE) <= 0 && mass.high.compareTo(ONE) >= 0)) {
    refuse("not-a-density", `the total mass is ${mass}, not 1; pass { normalize: true } to divide by it`);
  }
  return build(pieces, x, {
    name: "density",
    params: { f: pieces.map((p) => format(p.ast)).join(" | ") },
    constant,
    checks,
  });
}
/** f ≥ 0 on [lo, hi] by interval evaluation on ever smaller pieces: true, "negative", or null when undecided. */
function nonnegativeOn(ast, x, lo, hi, budget = 4000) {
  const stack = [[lo, hi]];
  let used = 0;
  while (stack.length) {
    const [a, b] = stack.pop();
    if (++used > budget) return null;
    let v;
    try {
      v = encloseAst(ast, { [x]: ivOf(a, b) }, 64);
    } catch {
      v = null;
    }
    if (v && v.lo.numerator >= 0n) continue;
    if (v && v.hi.numerator < 0n) return "negative";
    const m = a.add(b).divide(Q(2));
    if (b.subtract(a).compareTo(new Rational(1n, 1n << 40n)) < 0) return null;
    stack.push([a, m], [m, b]);
  }
  return true;
}

/** ∫_a^b f for a formula, certified (an end may be null for ±∞). */
function integralOf(ast, x, a, b, digits = 30) {
  const F0 = fn(ast, { variable: x });
  const r = a === null || b === null ? F0.improperIntegral(a ?? "-inf", b ?? "inf", { digits }).exact() : F0.integral(a, b, { digits }).exact();
  if (typeof r === "number") refuse("diverges", `the integral of ${format(ast)} diverges`);
  return tidyIv(r);
}
const times = (a, b) => ({ type: "op", op: "*", args: [a, b] });
const X = (x) => ({ type: "var", name: x });

/** The Continuous distribution for checked pieces. */
function build(pieces, x, extra) {
  const allPoly = pieces.every((p) => p.poly);
  // the CDF at the start of each piece
  const starts = [];
  let acc = ZERO;
  for (const p of pieces) {
    starts.push(acc);
    acc = add(acc, p.poly ? pIntegral(p.poly, p.a, p.b) : integralOf(p.ast, x, p.a, p.b));
  }
  const pieceAt = (q) => {
    for (let i = 0; i < pieces.length; i++) {
      const p = pieces[i];
      if ((p.a === null || q.compareTo(p.a) >= 0) && (p.b === null || q.compareTo(p.b) <= 0)) return i;
    }
    return -1;
  };
  const moment = (k) => {
    try {
      let s = ZERO;
      for (const p of pieces) {
        if (p.poly) s = add(s, pIntegral(pMul(p.poly, [...Array(k).fill(ZERO), ONE]), p.a, p.b));
        else {
          let t = p.ast;
          for (let i = 0; i < k; i++) t = times(X(x), t);
          s = add(s, integralOf(t, x, p.a, p.b, 24));
        }
      }
      return s;
    } catch {
      return null;
    }
  };
  const m1 = moment(1);
  const m2 = m1 === null ? null : moment(2);
  const mean = m1;
  const variance = m2 === null ? null : sub(m2, mul(m1, m1));
  const low = pieces[0].a;
  const high = pieces.at(-1).b;
  const fl = pieces.map((p) => (p.poly ? pFloat(p.poly) : ((F0) => (t) => F0({ [x]: t }))(compileFloat(p.ast))));
  const pdfAt = (q, bits) => {
    const i = pieceAt(q);
    if (i < 0) return ZERO;
    const p = pieces[i];
    return p.poly ? pEval(p.poly, q) : tidyIv(encloseAst(p.ast, { [x]: ivOf(q) }, bits).toRationalInterval());
  };
  const cdfAt = (q, bits) => {
    const i = pieceAt(q);
    const p = pieces[i];
    if (p.poly) return add(starts[i], pIntegral(p.poly, p.a, q));
    if (p.a !== null && q.equals(p.a)) return starts[i];
    return add(starts[i], integralOf(p.ast, x, p.a, q, Math.max(12, Math.ceil(bits / 3.3))));
  };
  // floating point: exact polynomials, or the certified value at low precision
  const fcdf = allPoly
    ? (t) => {
        if (low !== null && t <= low.toNumber()) return 0;
        if (high !== null && t >= high.toNumber()) return 1;
        const q = exactOf(t);
        const v = cdfAt(q, 53);
        return lowOf(v).add(highOf(v)).toNumber() / 2;
      }
    : (t) => {
        if (low !== null && t <= low.toNumber()) return 0;
        if (high !== null && t >= high.toNumber()) return 1;
        const v = cdfAt(exactOf(t), 40);
        return Math.min(1, Math.max(0, lowOf(v).add(highOf(v)).toNumber() / 2));
      };
  const [lf, hf] = [low === null ? -Infinity : low.toNumber(), high === null ? Infinity : high.toNumber()];
  const start = mean === null ? 0 : lowOf(mean).toNumber();
  return new Continuous({
    ...extra,
    low,
    high,
    mean,
    variance,
    _pieces: allPoly ? pieces.map((p) => ({ a: p.a, b: p.b, poly: p.poly })) : null,
    _asts: pieces.map((p) => ({ a: p.a, b: p.b, ast: p.ast })),
    pieces: pieces.map((p) => ({ from: p.a, to: p.b, density: p.poly ? pText(p.poly) : format(p.ast) })),
    _start: start,
    _pdf: (q, bits) => (q instanceof RationalInterval ? hullPdf(pieces, x, q, bits) : pdfAt(q, bits)),
    _edge: (where, bits) => pdfAt(where < 0 ? low : high, bits),
    _cdf: (q, bits) => cdfAt(q, bits),
    _fpdf: (t) => {
      if (t < lf || t > hf) return 0;
      const i = pieceAt(exactOf(t));
      return i < 0 ? 0 : fl[i](t);
    },
    _fcdf: fcdf,
    _draw: (rng) => F.quantileBy(fcdf, u01(rng), lf, hf, start),
  });
}
/** the density over an interval inside one piece, by interval arithmetic */
function hullPdf(pieces, x, iv, bits) {
  const p = pieces.find((pc) => (pc.a === null || iv.low.compareTo(pc.a) >= 0) && (pc.b === null || iv.high.compareTo(pc.b) <= 0));
  if (!p) refuse("unsupported", "the interval crosses a break between pieces of the density");
  return tidyIv(encloseAst(p.ast, { [x]: ivOf(iv.low, iv.high) }, bits).toRationalInterval());
}

// ---- transformations ------------------------------------------------------------------------

/** The pieces of a distribution as exact polynomials, when it has them. */
function polyPieces(dist) {
  if (dist._pieces) return dist._pieces;
  if (dist.name === "uniform") return [{ a: dist.low, b: dist.high, poly: [dist.high.subtract(dist.low).reciprocal()] }];
  return null;
}

/**
 * Y = g(X) for g proven monotone on the support of X (by the sign of g′ over
 * the whole support). The support must be bounded unless g is affine, and the
 * images of its ends must be exact. P(Y ≤ y) = P(X ≤ g⁻¹(y)) for increasing g
 * and P(X ≥ g⁻¹(y)) for decreasing g, with g⁻¹(y) enclosed by bisection;
 * the density is f(g⁻¹(y))/|g′(g⁻¹(y))|. The mean and variance are exact
 * when X has a piecewise polynomial density and g is a polynomial.
 */
export function transform(dist, g, options = {}) {
  const x = options.variable ?? "x";
  const ast = astOf(g);
  const vars = [...freeVariables(ast)];
  if (vars.some((v) => v !== x)) refuse("bad-argument", `g may only use the variable ${x}`);
  const dAst = simplify(differentiate(ast, x));
  const dConst = dAst.type === "num" ? dAst.value : null;
  if ((dist.low === null || dist.high === null) && dConst === null) refuse("unsupported", "on an unbounded support only an affine g is supported in this version");
  let sign;
  if (dConst !== null) sign = dConst.numerator > 0n ? 1 : dConst.numerator < 0n ? -1 : 0;
  else {
    let d;
    try {
      d = encloseAst(dAst, { [x]: ivOf(dist.low, dist.high) }, 96);
    } catch (e) {
      refuse("not-monotone", `g′ could not be enclosed on the support: ${e.message}`);
    }
    sign = d.lo.numerator >= 0n && !d.isPoint ? 1 : d.hi.numerator <= 0n && !d.isPoint ? -1 : 0;
    if (sign === 0 || (d.lo.numerator === 0n && d.hi.numerator === 0n)) refuse("not-monotone", `g = ${format(ast)} is not proven monotone on the support (g′ in [${d.lo}, ${d.hi}]); non-monotone g is not supported yet`);
  }
  if (sign === 0) refuse("not-monotone", "a constant g gives no continuous distribution");
  const gAt = (q) => {
    const v = encloseAst(ast, { [x]: ivOf(q) }, 96);
    if (!v.isPoint) refuse("unsupported", `g(${q}) is not an exact number, so the support of g(X) has no exact end`);
    return v.lo;
  };
  const ga = dist.low === null ? null : gAt(dist.low);
  const gb = dist.high === null ? null : gAt(dist.high);
  const [low, high] = sign > 0 ? [ga, gb] : [gb, ga];
  const gFloat = compileFloat(ast);
  const gf = (t) => gFloat({ [x]: t });
  const dFloat = compileFloat(dAst);
  // g⁻¹(y): the crossing of g (or −g, when decreasing) with y
  const inverse = (y, bits) => {
    const G = sign > 0 ? (t, b) => tidyIv(encloseAst(ast, { [x]: ivOf(t) }, b).toRationalInterval()) : (t, b) => tidyIv(encloseAst({ type: "op", op: "neg", args: [ast] }, { [x]: ivOf(t) }, b).toRationalInterval());
    const target = sign > 0 ? y : y.negate();
    // affine: exact
    if (dConst !== null) return y.subtract(gAtAffine(ast, x)).divide(dConst);
    // a simple rational with g(r) = y exactly (x² = 1/4 at 1/2), else bisection
    const r = simpleNear(fInv(y.toNumber()));
    if (r && (dist.low === null || r.compareTo(dist.low) >= 0) && (dist.high === null || r.compareTo(dist.high) <= 0)) {
      try {
        const v = encloseAst(ast, { [x]: ivOf(r) }, 64);
        if (v.isPoint && v.lo.equals(y)) return r;
      } catch {}
    }
    return crossing(G, target, { guess: fInv(y.toNumber()), lo: dist.low, hi: dist.high, bits });
  };
  const fInv = (yv) => {
    let [a, b] = [dist.low.toNumber(), dist.high.toNumber()];
    for (let i = 0; i < 80; i++) {
      const m = (a + b) / 2;
      if ((gf(m) - yv) * sign < 0) a = m;
      else b = m;
    }
    return (a + b) / 2;
  };
  // exact moments: E[g(X)^k] with g a polynomial and X piecewise polynomial
  const gp = polynomialOf(ast, x);
  const pp = polyPieces(dist);
  let mean = null;
  let variance = null;
  if (gp && pp) {
    const E = (k) => pp.reduce((s, p) => s.add(pIntegral(pMul(p.poly, pPow(gp, k)), p.a, p.b)), ZERO);
    mean = E(1);
    variance = E(2).subtract(mean.multiply(mean));
  } else if (dConst !== null) {
    // affine g = a x + b: E[g(X)] = a μ + b, Var = a² σ²
    const b0 = gAtAffine(ast, x);
    if (dist.mean !== null) mean = add(mul(dist.mean, dConst), b0);
    if (dist.variance !== null) variance = mul(dist.variance, dConst.multiply(dConst));
  }
  const fl = low === null ? -Infinity : low.toNumber();
  const fh = high === null ? Infinity : high.toNumber();
  const fcdf = (t) => {
    if (t <= fl) return 0;
    if (t >= fh) return 1;
    const xi = dConst !== null ? (t - gf(0)) / dConst.toNumber() : fInv(t);
    return sign > 0 ? dist.fast.cdf(xi) : 1 - dist.fast.cdf(xi);
  };
  return new Continuous({
    name: "transform",
    params: { of: String(dist), g: format(ast) },
    low,
    high,
    mean,
    variance,
    _start: (fl + fh) / 2,
    _pdf(y, bits) {
      const h = inverse(y, bits);
      const fx = isQ(h) ? dist.pdf(h, { bits }) : dist._pdf(h, bits);
      const d = tidyIv(encloseAst(dAst, { [x]: isQ(h) ? ivOf(h) : ivOf(h.low, h.high) }, bits).toRationalInterval());
      const ad = isQ(d) ? d.abs() : d.low.numerator >= 0n ? d : interval(d.high.negate(), d.low.negate());
      return div(fx, ad);
    },
    _cdf: (y, bits) => {
      const h = inverse(y, bits + 8);
      return sign > 0 ? dist.cdf(h, { bits }) : dist.sf(h, { bits });
    },
    _fpdf: (t) => {
      if (t <= fl || t >= fh) return 0;
      const xi = dConst !== null ? (t - gf(0)) / dConst.toNumber() : fInv(t);
      return dist.fast.pdf(xi) / Math.abs(dFloat({ [x]: xi }));
    },
    _fcdf: fcdf,
    _draw: (rng) => gf(dist.sample(rng)),
  });
}
/** the continued-fraction convergent of t with denominator at most 10⁶ that is within 1e-9, or null */
function simpleNear(t) {
  if (!Number.isFinite(t)) return null;
  let [h0, h1, k0, k1] = [0n, 1n, 1n, 0n];
  let v = t;
  for (let i = 0; i < 40; i++) {
    const a = Math.floor(v);
    const A = BigInt(a);
    [h0, h1] = [h1, A * h1 + h0];
    [k0, k1] = [k1, A * k1 + k0];
    if (k1 > 1000000n) return null;
    if (Math.abs(Number(h1) / Number(k1) - t) < 1e-9) return new Rational(h1, k1);
    const f = v - a;
    if (f < 1e-12) return new Rational(h1, k1);
    v = 1 / f;
  }
  return null;
}
function gAtAffine(ast, x) {
  return encloseAst(ast, { [x]: ivOf(ZERO) }, 64).lo;
}

// ---- order statistics -----------------------------------------------------------------------

const binom = (n, k) => {
  let r = 1n;
  for (let i = 1n; i <= BigInt(k); i++) r = (r * (BigInt(n) - i + 1n)) / i;
  return r;
};

/**
 * The k-th smallest of n independent draws from a distribution, as a
 * distribution. P(X₍ₖ₎ ≤ x) = Σ_{j ≥ k} C(n, j) F^j (1 − F)^(n−j), which is
 * increasing in F = F(x), so it is evaluated at the two ends of F's
 * enclosure; the density is n!/((k−1)!(n−k)!) F^(k−1) (1 − F)^(n−k) f. For a
 * piecewise polynomial density (uniform included) the result is again one,
 * with exact moments: the maximum of n uniforms is Beta(n, 1).
 */
export function orderStatistic(dist, k, n) {
  if (!Number.isInteger(n) || n < 1 || !Number.isInteger(k) || k < 1 || k > n) refuse("bad-argument", "orderStatistic(k, n) needs whole numbers 1 ≤ k ≤ n");
  const c = new Rational(binom(n, k) * BigInt(k));
  const B = (p) => {
    let s = ZERO;
    for (let j = k; j <= n; j++) s = s.add(new Rational(binom(n, j)).multiply(p.pow(j)).multiply(ONE.subtract(p).pow(n - j)));
    return s;
  };
  const Bf = (p) => {
    let s = 0;
    for (let j = k; j <= n; j++) s += Number(binom(n, j)) * p ** j * (1 - p) ** (n - j);
    return s;
  };
  const pp = polyPieces(dist);
  if (pp) {
    // the density c·F^(k−1)·(1 − F)^(n−k)·f, piece by piece, all exact
    let below = ZERO;
    const pieces = pp.map((p) => {
      const Fi = pAdd([below.subtract(pEval(pInt(p.poly), p.a))], pInt(p.poly));
      below = below.add(pIntegral(p.poly, p.a, p.b));
      const poly = pScale(pMul(pMul(pPow(Fi, k - 1), pPow(pAdd([ONE], pScale(Fi, Q(-1))), n - k)), p.poly), c);
      return { a: p.a, b: p.b, poly, ast: parse(pText(poly)) };
    });
    return build(pieces, "x", { name: "orderStatistic", params: { k, n, of: String(dist) }, constant: ONE, checks: { nonnegative: "by construction", mass: "exact" } });
  }
  const g = (p) => p.pow(k - 1).multiply(ONE.subtract(p).pow(n - k));
  // the range of F^(k−1)(1 − F)^(n−k) over an enclosure of F: it rises to its peak at (k−1)/(n−1), then falls
  const gRange = (Fv) => {
    if (isQ(Fv)) return g(Fv);
    const [a, b] = [g(Fv.low), g(Fv.high)];
    const lo = a.compareTo(b) <= 0 ? a : b;
    let hi = a.compareTo(b) >= 0 ? a : b;
    if (n > 1) {
      const peak = new Rational(BigInt(k - 1), BigInt(n - 1));
      if (peak.compareTo(Fv.low) > 0 && peak.compareTo(Fv.high) < 0) hi = g(peak);
    }
    return interval(lo, hi);
  };
  return new Continuous({
    name: "orderStatistic",
    params: { k, n, of: String(dist) },
    low: dist.low,
    high: dist.high,
    mean: null,
    variance: null,
    _start: dist._start ?? 0,
    _pdf: (q, bits) => mul(c, mul(gRange(dist.cdf(q, { bits: bits + 8 })), dist.pdf(q, { bits: bits + 8 }))),
    _cdf: (q, bits) => {
      const Fv = dist.cdf(q, { bits: bits + 8 });
      return isQ(Fv) ? B(Fv) : interval(B(Fv.low), B(Fv.high));
    },
    _fpdf: (t) => {
      const p = dist.fast.cdf(t);
      return Number(c.numerator) * p ** (k - 1) * (1 - p) ** (n - k) * dist.fast.pdf(t);
    },
    _fcdf: (t) => Bf(dist.fast.cdf(t)),
    _draw(rng) {
      const xs = dist.samples(n, rng).sort((a, b) => a - b);
      return xs[k - 1];
    },
  });
}

// every distribution gets the two constructions as methods
Continuous.prototype.transform = function (g, options) {
  return transform(this, g, options);
};
Continuous.prototype.orderStatistic = function (k, n) {
  return orderStatistic(this, k, n);
};
