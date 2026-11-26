// Numbers for the statistics library: reading what people pass in as exact
// rationals or rational intervals, a little interval arithmetic, and the
// three kinds of answer (exact, certified, simulated) and how they are shown.
import { Rational, RationalInterval, parseNumber } from "@ratmath/core";
import { enclose } from "@ratmath/analytic";

export class StatsError extends Error {
  constructor(code, message) {
    super(message);
    this.name = "StatsError";
    this.code = code;
  }
}
export const refuse = (code, message) => {
  throw new StatsError(code, message);
};

export const Q = (n, d = 1n) => new Rational(BigInt(n), BigInt(d));
export const ZERO = Q(0);
export const ONE = Q(1);
export const TWO = Q(2);
export const HALF = Q(1, 2);

export const isQ = (v) => v instanceof Rational;
export const isIv = (v) => v instanceof RationalInterval;

/** A JS number as the exact rational its shortest decimal form names: 0.1 is 1/10. */
function numberToQ(v, what) {
  if (!Number.isFinite(v)) refuse("not-a-number", `${what} must be a finite number, got ${v}`);
  if (Number.isInteger(v) && Math.abs(v) < 1e21) return new Rational(BigInt(v));
  const text = String(v);
  const m = /^(-?)(\d+)(?:\.(\d+))?(?:e([+-]\d+))?$/.exec(text);
  /* c8 ignore next */
  if (!m) refuse("not-a-number", `cannot read ${text} as a number`);
  const digits = BigInt(m[2] + (m[3] ?? ""));
  const exp = Number(m[4] ?? 0) - (m[3]?.length ?? 0);
  const q = exp >= 0 ? new Rational(digits * 10n ** BigInt(exp)) : new Rational(digits, 10n ** BigInt(-exp));
  return m[1] ? q.negate() : q;
}

/** Read an exact rational: a Rational, BigInt, JS number (by its decimal form) or text such as "1/3" or "0.25". */
export function toQ(v, what = "a value") {
  if (v instanceof Rational) return v;
  if (v instanceof RationalInterval) {
    if (v.low.equals(v.high)) return v.low;
    refuse("needs-exact", `${what} must be an exact number here, not an interval`);
  }
  if (typeof v === "bigint") return new Rational(v);
  if (typeof v === "number") return numberToQ(v, what);
  if (typeof v === "string") {
    const t = v.trim();
    const sci = /^(-?\d+(?:\.\d+)?)e([+-]?\d+)$/i.exec(t);
    try {
      if (sci) {
        const e = Number(sci[2]);
        const scale = e >= 0 ? new Rational(10n ** BigInt(e)) : new Rational(1n, 10n ** BigInt(-e));
        return Rational.from(sci[1]).multiply(scale);
      }
      return Rational.from(t);
    } catch {
      refuse("not-a-number", `cannot read "${v}" as an exact number for ${what}`);
    }
  }
  if (v && typeof v === "object" && typeof v.numerator === "bigint" && typeof v.denominator === "bigint") return new Rational(v.numerator, v.denominator);
  if (v && typeof v === "object" && typeof v.value === "bigint") return new Rational(v.value);
  return refuse("not-a-number", `${what} must be a number, got ${describe(v)}`);
}
const describe = (v) => (v === null || v === undefined ? String(v) : typeof v === "object" ? (v.constructor?.name ?? "object") : `${typeof v} ${String(v)}`);

/**
 * Read a data value: an exact rational, or an interval for a measurement with
 * stated uncertainty (a RationalInterval, a pair [low, high], or Core's
 * notations "1.2:1.4" and "1.23[+-0.05]").
 */
export function toValue(v, what = "a data value") {
  if (v instanceof RationalInterval) return v.low.equals(v.high) ? v.low : v;
  if (Array.isArray(v) && v.length === 2) return interval(toQ(v[0], what), toQ(v[1], what));
  if (typeof v === "string" && /[:[]/.test(v)) {
    try {
      const parsed = parseNumber(v.trim());
      if (parsed instanceof RationalInterval) return parsed.low.equals(parsed.high) ? parsed.low : parsed;
      return toQ(parsed, what);
    } catch {
      refuse("not-a-number", `cannot read "${v}" as a number or an interval for ${what}`);
    }
  }
  return toQ(v, what);
}

export const interval = (lo, hi) => (lo.compareTo(hi) <= 0 ? new RationalInterval(lo, hi) : new RationalInterval(hi, lo));
export const lowOf = (v) => (v instanceof RationalInterval ? v.low : v);
export const highOf = (v) => (v instanceof RationalInterval ? v.high : v);
const asIv = (v) => (v instanceof RationalInterval ? v : new RationalInterval(v, v));
const tidy = (iv) => (iv.low.equals(iv.high) ? iv.low : iv);

// arithmetic on values that are rationals or intervals
export const add = (a, b) => (isQ(a) && isQ(b) ? a.add(b) : tidy(asIv(a).add(asIv(b))));
export const sub = (a, b) => (isQ(a) && isQ(b) ? a.subtract(b) : tidy(asIv(a).subtract(asIv(b))));
export const mul = (a, b) => (isQ(a) && isQ(b) ? a.multiply(b) : tidy(asIv(a).multiply(asIv(b))));
export const div = (a, b) => {
  if (isQ(a) && isQ(b)) {
    if (b.numerator === 0n) refuse("division-by-zero", "division by zero");
    return a.divide(b);
  }
  const d = asIv(b);
  if (d.low.numerator <= 0n && d.high.numerator >= 0n) refuse("division-by-zero", `cannot divide by an interval that contains 0 (${d})`);
  return tidy(asIv(a).divide(d));
};
export const neg = (a) => (isQ(a) ? a.negate() : new RationalInterval(a.high.negate(), a.low.negate()));
/** The square of a value: for an interval, the exact range of x². */
export const square = (a) => {
  if (isQ(a)) return a.multiply(a);
  const l = a.low.multiply(a.low);
  const h = a.high.multiply(a.high);
  if (a.low.numerator <= 0n && a.high.numerator >= 0n) return interval(ZERO, max(l, h));
  return interval(min(l, h), max(l, h));
};
export const min = (a, b) => (a.compareTo(b) <= 0 ? a : b);
export const max = (a, b) => (a.compareTo(b) >= 0 ? a : b);
export const cmp = (a, b) => a.compareTo(b);

const isqrt = (n) => {
  if (n < 2n) return n;
  let x = 1n << BigInt((n.toString(2).length + 1) >> 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
};
/** √q: the exact rational when q is a perfect square, otherwise a certified enclosure. */
export function sqrt(v, bits = 96) {
  if (!isQ(v) && !(v instanceof RationalInterval)) v = toValue(v, "the number");
  if (isQ(v)) {
    if (v.numerator < 0n) refuse("domain", `the square root of a negative number (${v})`);
    const n = isqrt(v.numerator);
    const d = isqrt(v.denominator);
    if (n * n === v.numerator && d * d === v.denominator) return new Rational(n, d);
    return enclose.sqrt(v, bits);
  }
  const lo = v.low.numerator < 0n ? ZERO : v.low;
  if (v.high.numerator < 0n) refuse("domain", `the square root of a negative interval (${v})`);
  return interval(lowOf(sqrt(lo, bits)), highOf(sqrt(v.high, bits)));
}

/** The midpoint of a value as a JS number. */
export const toNumber = (v) => (isQ(v) ? v.toNumber() : v instanceof RationalInterval ? v.low.add(v.high).divide(TWO).toNumber() : Number(v));

const pow10 = (n) => 10n ** BigInt(n);
/** q rounded to `digits` decimal places, as text (half away from zero). */
function fixed(q, digits) {
  const scaled = q.multiply(new Rational(pow10(digits)));
  const twice = scaled.multiply(TWO);
  const neg0 = twice.numerator < 0n;
  const a = neg0 ? twice.negate() : twice;
  const r = (a.numerator / a.denominator + 1n) / 2n;
  let s = r.toString().padStart(digits + 1, "0");
  if (digits > 0) s = `${s.slice(0, -digits)}.${s.slice(-digits)}`;
  return neg0 && r !== 0n ? `-${s}` : s;
}
/**
 * A value as decimal text. A rational is rounded to `digits` places. An
 * interval is shown by the digits its two ends share, followed by "…", or as
 * [low, high] when they share none at that many places.
 */
export function show(v, digits = 6) {
  if (typeof v === "number") return String(Number(v.toPrecision(digits + 1)));
  if (isQ(v)) return v.denominator === 1n ? v.toString() : fixed(v, digits);
  if (!(v instanceof RationalInterval)) return String(v);
  if (v.low.equals(v.high)) return show(v.low, digits);
  for (let d = digits; d >= 0; d--) {
    const a = fixed(v.low, d);
    const b = fixed(v.high, d);
    if (a === b) return d === digits ? a : `${a}…`;
  }
  return `[${fixed(v.low, digits)}, ${fixed(v.high, digits)}]`;
}

/** Is the value known to be below, above, or not separable from the level? */
export function compare(value, level) {
  const l = toQ(level, "the level");
  if (highOf(value).compareTo(l) < 0) return "below";
  if (lowOf(value).compareTo(l) > 0) return "above";
  if (isQ(value)) return "equal";
  return "too close to call";
}

/** k as a non-negative whole number (a JS number), or a refusal. */
export function whole(v, what, { min: least = 0 } = {}) {
  const q = typeof v === "number" || typeof v === "bigint" ? v : toQ(v, what);
  const n = typeof q === "object" ? (q.denominator === 1n ? Number(q.numerator) : NaN) : Number(q);
  if (!Number.isSafeInteger(n) || n < least) refuse("bad-argument", `${what} must be a whole number${least > 0 ? ` of at least ${least}` : least === 0 ? " that is not negative" : ""}, got ${String(v)}`);
  return n;
}
/** p as an exact probability in [0, 1]. */
export function probability(v, what = "a probability") {
  const q = toQ(v, what);
  if (q.numerator < 0n || q.compareTo(ONE) > 0) refuse("bad-argument", `${what} must lie between 0 and 1, got ${q}`);
  return q;
}
