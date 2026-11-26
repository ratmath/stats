// A certified quantile: the point where an increasing function with
// certified values crosses a level, enclosed between two rationals.
import { Rational, RationalInterval } from "@ratmath/core";
import { refuse, lowOf, highOf, Q, TWO } from "./num.js";

const pow2 = (k) => (k >= 0 ? new Rational(1n << BigInt(k)) : new Rational(1n, 1n << BigInt(-k)));
/** q rounded to a multiple of 2^(−bits). */
const dyadic = (q, bits) => {
  if (bits <= 0) return new Rational(q.floor());
  const scale = 1n << BigInt(bits);
  return new Rational(q.multiply(new Rational(scale)).floor(), scale);
};
const fromFloat = (x) => {
  // the exact value of a double
  if (x === 0) return Q(0);
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  const hi = buf.getUint32(0);
  const lo = buf.getUint32(4);
  const exp = (hi >>> 20) & 0x7ff;
  const mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  const m = exp === 0 ? mant : mant | (1n << 52n);
  const e = (exp === 0 ? 1 : exp) - 1075;
  const v = e >= 0 ? new Rational(m << BigInt(e)) : new Rational(m, 1n << BigInt(-e));
  return hi >>> 31 ? v.negate() : v;
};

/**
 * Enclose the x with F(x) = p, for F increasing on (lo, hi).
 *   F(x, bits)   a certified enclosure of F at the rational x
 *   p            the level, a rational strictly between F's limits
 *   guess        a float near the answer
 *   lo, hi       the ends of the support as rationals, or null for infinite
 * The result is an interval of width at most 2^(−bits) relative to the size
 * of the answer.
 */
export function crossing(F, p, { guess, lo = null, hi = null, bits = 64 }) {
  // The sign of F(x) − p, with the midpoint of the enclosure that settled it.
  // The precision starts where the last call needed it, since the points
  // close in on the crossing and need more each time.
  let level = 32;
  let last = null;
  const side = (x) => {
    for (;;) {
      const v = F(x, level);
      last = lowOf(v).add(highOf(v)).divide(TWO).subtract(p);
      if (highOf(v).compareTo(p) < 0) return -1;
      if (lowOf(v).compareTo(p) > 0) return 1;
      if (level >= 2 * bits + 48) return 0;
      level = Math.min(2 * bits + 48, level + 32);
    }
  };
  if (!Number.isFinite(guess)) guess = lo !== null ? lo.toNumber() + 1 : hi !== null ? hi.toNumber() - 1 : 0;
  const g = fromFloat(guess);
  const size = Math.max(Math.abs(guess), 2 ** -40);
  const scaleBits = Math.ceil(Math.log2(size));
  const tol = pow2(scaleBits - bits);
  const grid = bits - scaleBits + 6;
  // a first bracket around the guess, widened until it holds the crossing;
  // next to an end of the support the candidates close in on that end instead
  let a = null;
  let b = null;
  let fa = null;
  let fb = null;
  let step = pow2(scaleBits - 24);
  let centre = g;
  if (lo !== null && centre.compareTo(lo) <= 0) centre = lo.add(step);
  if (hi !== null && centre.compareTo(hi) >= 0) centre = hi.subtract(step);
  let ka = 1;
  let kb = 1;
  for (let i = 0; i < 600 && (a === null || b === null); i++) {
    if (a === null) {
      let x = centre.subtract(step);
      if (lo !== null && x.compareTo(lo) <= 0) x = lo.add(centre.subtract(lo).divide(pow2(ka++)));
      if (side(x) < 0) {
        a = x;
        fa = last;
      }
    }
    if (b === null) {
      let x = centre.add(step);
      if (hi !== null && x.compareTo(hi) >= 0) x = hi.subtract(hi.subtract(centre).divide(pow2(kb++)));
      if (side(x) > 0) {
        b = x;
        fb = last;
      }
    }
    step = step.multiply(Q(16));
  }
  if (a === null || b === null) refuse("no-bracket", "could not bracket the quantile; the level may be too close to 0 or 1 for this precision");

  let slow = 0;
  for (let i = 0; i < 400; i++) {
    const width = b.subtract(a);
    if (width.compareTo(tol) <= 0) break;
    // a secant step, pulled a little toward the far end so that both ends move; bisection if progress stalls
    let c;
    const denom = fb.subtract(fa);
    if (slow < 2 && denom.numerator > 0n) {
      c = a.subtract(fa.multiply(width).divide(denom));
      const nudge = tol.divide(Q(4));
      c = c.subtract(a).compareTo(b.subtract(c)) < 0 ? c.add(nudge) : c.subtract(nudge);
      if (c.compareTo(a) <= 0 || c.compareTo(b) >= 0) c = a.add(b).divide(TWO);
    } else {
      c = a.add(b).divide(TWO);
      slow = 0;
    }
    c = dyadic(c, grid);
    if (c.compareTo(a) <= 0 || c.compareTo(b) >= 0) {
      c = a.add(b).divide(TWO);
    }
    const s = side(c);
    const before = width;
    if (s < 0) {
      a = c;
      fa = last;
    } else if (s > 0) {
      b = c;
      fb = last;
    } else {
      // F(c) is p as far as can be told: step half a tolerance to each side
      const h = tol.divide(TWO);
      const left = c.subtract(h);
      const right = c.add(h);
      if (left.compareTo(a) > 0 && side(left) < 0) a = left;
      if (right.compareTo(b) < 0 && side(right) > 0) b = right;
      if (b.subtract(a).compareTo(before) >= 0) refuse("too-close", "the quantile could not be separated at this precision");
      break;
    }
    if (b.subtract(a).multiply(TWO).compareTo(before) > 0) slow++;
    else slow = 0;
  }
  return new RationalInterval(a, b);
}
