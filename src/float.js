// Floating-point special functions for the fast side of the library:
// previews, first guesses for certified root finding, and simulation. None
// of these is certified; the certified versions are in @ratmath/analytic.
import { erfFloat } from "@ratmath/analytic";

const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
/** log Γ(x) for x > 0. */
export function lgamma(x) {
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
  const y = x - 1;
  let a = LANCZOS[0];
  const t = y + 7.5;
  for (let i = 1; i < 9; i++) a += LANCZOS[i] / (y + i);
  return 0.5 * Math.log(2 * Math.PI) + (y + 0.5) * Math.log(t) - t + Math.log(a);
}
export const erf = erfFloat;
/** Φ(z), the standard normal cumulative probability, accurate in both tails. */
export function phi(z) {
  if (z < -38) return 0;
  if (z > 38) return 1;
  const a = Math.abs(z) / Math.SQRT2;
  const tail = a < 3 ? 0.5 * (1 - erfFloat(a)) : 0.5 * erfcBig(a);
  return z < 0 ? tail : 1 - tail;
}
/** erfc(a) for a ≥ 3 by its continued fraction. */
function erfcBig(a) {
  let f = 0;
  for (let k = 60; k >= 1; k--) f = k / 2 / (a + f);
  return Math.exp(-a * a) / (Math.sqrt(Math.PI) * (a + f));
}
/** The regularized incomplete gamma function P(a, x). */
export function gammaP(a, x) {
  if (x <= 0) return 0;
  const lead = Math.exp(a * Math.log(x) - x - lgamma(a));
  if (x < a + 1) {
    let term = 1 / a;
    let s = term;
    for (let n = 1; n < 100000; n++) {
      term *= x / (a + n);
      s += term;
      if (term < s * 1e-17) break;
    }
    return Math.min(1, lead * s);
  }
  // Lentz's method for the continued fraction of Q(a, x)
  const tiny = 1e-300;
  let b = x + 1 - a;
  let c = 1 / tiny;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 100000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return Math.max(0, 1 - lead * h);
}
/** The regularized incomplete beta function I_x(a, b). */
export function betaI(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  if (x > (a + 1) / (a + b + 2)) return 1 - betaI(1 - x, b, a);
  const lead = Math.exp(a * Math.log(x) + b * Math.log1p(-x) + lgamma(a + b) - lgamma(a) - lgamma(b)) / a;
  const tiny = 1e-300;
  let c = 1;
  let d = 1 - ((a + b) * x) / (a + 1);
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let h = d;
  for (let m = 1; m < 100000; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((a + m2 - 1) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (a + b + m) * x) / ((a + m2) * (a + m2 + 1));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return Math.min(1, lead * h);
}
const A = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239];
const Bc = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
const C = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
const D = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];
/** The standard normal quantile: Acklam's rational approximation, then one Halley step. */
export function normalQuantile(p) {
  if (!(p > 0)) return -Infinity;
  if (!(p < 1)) return Infinity;
  let x;
  if (p < 0.02425) {
    const q = Math.sqrt(-2 * Math.log(p));
    x = (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1);
  } else if (p > 1 - 0.02425) {
    const q = Math.sqrt(-2 * Math.log1p(-p));
    x = -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1);
  } else {
    const q = p - 0.5;
    const r = q * q;
    x = ((((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q) / (((((Bc[0] * r + Bc[1]) * r + Bc[2]) * r + Bc[3]) * r + Bc[4]) * r + 1);
  }
  const e = phi(x) - p;
  const u = e * Math.sqrt(2 * Math.PI) * Math.exp((x * x) / 2);
  return x - u / (1 + (x * u) / 2);
}
/**
 * The p-quantile of an increasing cumulative function by bracketing and
 * bisection, on [lo, hi] (either may be infinite).
 */
export function quantileBy(cdf, p, lo, hi, start = 0) {
  if (!(p > 0)) return lo;
  if (!(p < 1)) return hi;
  let a = Number.isFinite(lo) ? lo : Math.min(start, Number.isFinite(hi) ? hi : start) - 1;
  let b = Number.isFinite(hi) ? hi : Math.max(start, a) + 1;
  for (let step = 1; !Number.isFinite(lo) && cdf(a) > p && step < 1e300; step *= 2) a -= step;
  for (let step = 1; !Number.isFinite(hi) && cdf(b) < p && step < 1e300; step *= 2) b += step;
  for (let i = 0; i < 200; i++) {
    const m = (a + b) / 2;
    if (m === a || m === b) break;
    if (cdf(m) < p) a = m;
    else b = m;
  }
  return (a + b) / 2;
}
