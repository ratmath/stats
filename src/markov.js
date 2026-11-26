// Finite Markov chains, exactly. A chain is a square matrix of rational
// transition probabilities with rows summing to 1. Everything a course asks
// of one is a finite computation over the rationals: n-step distributions,
// communicating classes and their periods, stationary distributions,
// absorption probabilities, and expected hitting times. Paths are simulated
// from a seed.
import { Rational } from "@ratmath/core";
import { rng as makeRng } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, Q, whole, probability, show } from "./num.js";

const gcd = (a, b) => (b === 0 ? a : gcd(b, a % b));

/** Solve A·x = b exactly (A square, rationals); null when A is singular. */
function solve(A, b) {
  const n = A.length;
  if (n === 0) return [];
  const M = A.map((row, i) => [...row, ...(Array.isArray(b[i]) ? b[i] : [b[i]])]);
  const w = M[0].length;
  for (let c = 0; c < n; c++) {
    let p = c;
    while (p < n && M[p][c].numerator === 0n) p++;
    if (p === n) return null;
    [M[c], M[p]] = [M[p], M[c]];
    const inv = ONE.divide(M[c][c]);
    for (let j = c; j < w; j++) M[c][j] = M[c][j].multiply(inv);
    for (let r = 0; r < n; r++) {
      if (r === c || M[r][c].numerator === 0n) continue;
      const f = M[r][c];
      for (let j = c; j < w; j++) M[r][j] = M[r][j].subtract(f.multiply(M[c][j]));
    }
  }
  return M.map((row) => (w === n + 1 ? row[n] : row.slice(n)));
}

export class MarkovChain {
  constructor(P, states) {
    this.P = P;
    this.states = states;
    this.size = P.length;
    Object.freeze(this);
  }
  #at(s) {
    const i = this.states.findIndex((t) => t === s || String(t) === String(s));
    if (i < 0) refuse("bad-argument", `${s} is not a state of this chain`);
    return i;
  }
  #dist(d) {
    if (d === undefined) refuse("bad-argument", "a starting state or distribution is needed");
    if (Array.isArray(d)) {
      if (d.length !== this.size) refuse("bad-argument", `a distribution needs ${this.size} probabilities`);
      const v = d.map((p) => probability(p, "a probability"));
      if (!v.reduce((s, p) => s.add(p), ZERO).equals(ONE)) refuse("bad-argument", "a distribution must add up to 1");
      return v;
    }
    if (typeof d === "object") return this.states.map((s) => probability(d[s] ?? 0, "a probability"));
    const v = new Array(this.size).fill(ZERO);
    v[this.#at(d)] = ONE;
    return v;
  }
  /** The probability of going from state a to state b in one step. */
  p(a, b) {
    return this.P[this.#at(a)][this.#at(b)];
  }
  /** The n-step transition matrix Pⁿ. */
  power(n) {
    let e = whole(n, "n");
    let R = this.P.map((r, i) => r.map((_, j) => (i === j ? ONE : ZERO)));
    let B = this.P;
    const mul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((s, a, t) => (a.numerator === 0n ? s : s.add(a.multiply(Y[t][j]))), ZERO)));
    while (e > 0) {
      if (e & 1) R = mul(R, B);
      e >>= 1;
      if (e) B = mul(B, B);
    }
    return R;
  }
  /** The distribution after n steps from a state, a list of probabilities, or { state: p }. Keyed by state. */
  distribution(start, n = 1) {
    const v = this.#dist(start);
    const Pn = this.power(n);
    const out = {};
    this.states.forEach((s, j) => (out[s] = v.reduce((acc, p, i) => acc.add(p.multiply(Pn[i][j])), ZERO)));
    return out;
  }
  /**
   * The communicating classes, each with { states, closed, recurrent, period }.
   * In a finite chain a class is recurrent exactly when it is closed.
   */
  classes() {
    const n = this.size;
    const reach = this.P.map((row, i) => {
      const seen = new Set([i]);
      const q = [i];
      while (q.length) {
        const a = q.pop();
        this.P[a].forEach((p, b) => {
          if (p.numerator !== 0n && !seen.has(b)) {
            seen.add(b);
            q.push(b);
          }
        });
      }
      return seen;
    });
    const done = new Set();
    const out = [];
    for (let i = 0; i < n; i++) {
      if (done.has(i)) continue;
      const cls = [...reach[i]].filter((j) => reach[j].has(i)).sort((a, b) => a - b);
      cls.forEach((j) => done.add(j));
      const inside = new Set(cls);
      const closed = cls.every((a) => this.P[a].every((p, b) => p.numerator === 0n || inside.has(b)));
      // period: gcd of level differences along edges inside the class, from a breadth-first search
      const level = new Map([[cls[0], 0]]);
      const q = [cls[0]];
      let g = 0;
      for (let h = 0; h < q.length; h++) {
        const a = q[h];
        this.P[a].forEach((p, b) => {
          if (p.numerator === 0n || !inside.has(b)) return;
          if (!level.has(b)) {
            level.set(b, level.get(a) + 1);
            q.push(b);
          } else g = gcd(g, Math.abs(level.get(a) + 1 - level.get(b)));
        });
      }
      out.push(Object.freeze({ states: cls.map((j) => this.states[j]), closed, recurrent: closed, period: g || null }));
    }
    return out;
  }
  isIrreducible() {
    return this.classes().length === 1;
  }
  /** Irreducible with period 1: then Pⁿ converges to rows equal to the stationary distribution. */
  isErgodic() {
    const c = this.classes();
    return c.length === 1 && c[0].period === 1;
  }
  /**
   * The stationary distribution π = πP, Σπ = 1, exactly. A chain with more
   * than one closed class has one for each; they are all returned, keyed by
   * state, in `all`, and `unique` says whether there is just one.
   */
  stationary() {
    const closed = this.classes().filter((c) => c.closed);
    const all = closed.map((c) => {
      const idx = c.states.map((s) => this.#at(s));
      const m = idx.length;
      // (Pᵀ − I) restricted to the class, with the last equation replaced by Σπ = 1
      const A = idx.map((_, r) => idx.map((_, k) => this.P[idx[k]][idx[r]].subtract(r === k ? ONE : ZERO)));
      A[m - 1] = idx.map(() => ONE);
      const b = idx.map((_, r) => (r === m - 1 ? ONE : ZERO));
      const pi = solve(A, b);
      const out = {};
      this.states.forEach((s) => (out[s] = ZERO));
      idx.forEach((j, r) => (out[this.states[j]] = pi[r]));
      return out;
    });
    return Object.freeze({ unique: all.length === 1, distribution: all.length === 1 ? all[0] : null, all });
  }
  /**
   * Absorption: for each transient state, the probability of ending in each
   * closed class, and the expected number of steps before entering one.
   * Uses the fundamental matrix N = (I − Q)⁻¹, Q the transient-to-transient block.
   */
  absorption() {
    const cls = this.classes();
    const closed = cls.filter((c) => c.closed);
    const transient = cls.filter((c) => !c.closed).flatMap((c) => c.states).map((s) => this.#at(s)).sort((a, b) => a - b);
    const t = transient.length;
    const name = (c) => (c.states.length === 1 ? String(c.states[0]) : `{${c.states.join(", ")}}`);
    if (!t) return Object.freeze({ transient: [], targets: closed.map(name), probabilities: {}, steps: {}, fundamental: [] });
    const IQ = transient.map((a, r) => transient.map((b, k) => (r === k ? ONE : ZERO).subtract(this.P[a][b])));
    const R = transient.map((a) => closed.map((c) => c.states.map((s) => this.#at(s)).reduce((s, b) => s.add(this.P[a][b]), ZERO)));
    const I = transient.map((_, r) => transient.map((_, k) => (r === k ? ONE : ZERO)));
    const N = solve(IQ, I);
    const probabilities = {};
    const steps = {};
    transient.forEach((a, r) => {
      const row = {};
      closed.forEach((c, k) => (row[name(c)] = N[r].reduce((s, v, j) => s.add(v.multiply(R[j][k])), ZERO)));
      probabilities[this.states[a]] = row;
      steps[this.states[a]] = N[r].reduce((s, v) => s.add(v), ZERO);
    });
    return Object.freeze({ transient: transient.map((a) => this.states[a]), targets: closed.map(name), probabilities, steps, fundamental: N });
  }
  /**
   * Expected steps to reach a set of target states from each state (0 on the
   * targets; Infinity where the targets may never be reached), and the
   * probability of ever reaching them.
   */
  hitting(targets) {
    const T = new Set((Array.isArray(targets) ? targets : [targets]).map((s) => this.#at(s)));
    const n = this.size;
    // the states that can reach T at all
    const can = new Set(T);
    let grew = true;
    while (grew) {
      grew = false;
      for (let a = 0; a < n; a++) if (!can.has(a) && this.P[a].some((p, b) => p.numerator !== 0n && can.has(b))) {
        can.add(a);
        grew = true;
      }
    }
    const others = [...Array(n).keys()].filter((a) => !T.has(a) && can.has(a));
    const A = others.map((a, r) => others.map((b, k) => (r === k ? ONE : ZERO).subtract(this.P[a][b])));
    const prob = solve(A, others.map((a) => [...T].reduce((s, b) => s.add(this.P[a][b]), ZERO))) ?? [];
    const probability = {};
    const expected = {};
    this.states.forEach((s, a) => {
      probability[s] = T.has(a) ? ONE : can.has(a) ? prob[others.indexOf(a)] : ZERO;
    });
    // expected time, where the target is reached for sure
    const sure = others.filter((a) => probability[this.states[a]].equals(ONE));
    const B = sure.map((a, r) => sure.map((b, k) => (r === k ? ONE : ZERO).subtract(this.P[a][b])));
    const time = solve(B, sure.map(() => ONE)) ?? [];
    this.states.forEach((s, a) => {
      expected[s] = T.has(a) ? ZERO : sure.includes(a) ? time[sure.indexOf(a)] : Infinity;
    });
    return Object.freeze({ probability, expected });
  }
  /** The expected time to return to a state: 1/π(state) for an irreducible chain (Kac). */
  meanReturnTime(state) {
    const i = this.#at(state);
    const h = this.hitting([state]).expected;
    let s = ONE;
    for (let b = 0; b < this.size; b++) {
      const p = this.P[i][b];
      if (b === i || p.numerator === 0n) continue;
      const e = h[this.states[b]];
      if (e === Infinity) return Infinity;
      s = s.add(p.multiply(e));
    }
    return s;
  }
  /** A path of the given number of steps from a start state, drawn from a seed. */
  simulate(start, steps, { seed = 1 } = {}) {
    const k = whole(steps, "steps");
    const g = makeRng(seed).fork("Markov chain");
    let i = this.#at(start);
    const path = [this.states[i]];
    for (let t = 0; t < k; t++) {
      i = g.weightedIndex(this.P[i]);
      path.push(this.states[i]);
    }
    return path;
  }
  toString() {
    const w = Math.max(...this.states.map((s) => String(s).length), ...this.P.flat().map((p) => String(p).length));
    const pad = (v) => String(v).padStart(w);
    return [" ".repeat(w) + " " + this.states.map(pad).join(" "), ...this.P.map((r, i) => pad(this.states[i]) + " " + r.map(pad).join(" "))].join("\n");
  }
}

/**
 * A finite Markov chain from its transition matrix (rows: from, columns:
 * to), as a list of rows or an object { from: { to: p } }. Each row must add
 * up to 1 exactly.
 *   markovChain([[ "1/2", "1/2" ], [ "1/3", "2/3" ]], { states: ["sun", "rain"] })
 *   markovChain({ sun: { sun: "1/2", rain: "1/2" }, rain: { sun: "1/3", rain: "2/3" } })
 */
export function markovChain(matrix, { states } = {}) {
  let names;
  let rows;
  if (Array.isArray(matrix)) {
    rows = matrix;
    names = states ?? matrix.map((_, i) => i);
  } else if (matrix && typeof matrix === "object") {
    const keys = new Set(Object.keys(matrix));
    for (const r of Object.values(matrix)) for (const k of Object.keys(r)) keys.add(k);
    names = states ?? [...keys];
    rows = names.map((a) => names.map((b) => matrix[a]?.[b] ?? 0));
  } else refuse("bad-argument", "markovChain needs a matrix (a list of rows) or an object of rows");
  const n = rows.length;
  if (!n) refuse("no-data", "a chain needs at least one state");
  if (names.length !== n) refuse("bad-argument", `${n} rows but ${names.length} state names`);
  if (new Set(names.map(String)).size !== n) refuse("bad-argument", "state names must differ");
  const P = rows.map((r, i) => {
    if (!Array.isArray(r) || r.length !== n) refuse("bad-argument", "the matrix must be square");
    const row = r.map((p) => probability(p, "a transition probability"));
    const s = row.reduce((a, b) => a.add(b), ZERO);
    if (!s.equals(ONE)) refuse("bad-argument", `the row of ${names[i]} adds up to ${s}, not 1`);
    return Object.freeze(row);
  });
  return new MarkovChain(Object.freeze(P), Object.freeze(names.slice()));
}
