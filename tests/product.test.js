import { describe, expect, test } from "bun:test";
import * as S from "../index.js";

const str = (v) => String(v);
// The same distribution listed as a table, the way products used to be built.
const listed = (d) => new S.Finite(d.entries);
const rows = (d) => d.entries.map(([o, p]) => [JSON.stringify(o), str(p)]);
const loaded = () => S.finite([1, 2, 3, 4, 5, 6], ["1/12", "1/12", "1/12", "1/12", "1/3", "1/3"]);

const cases = {
  "fair dice": () => S.dice(3),
  "loaded dice": () => loaded().power(3),
  "mixed parts": () => S.coin("2/7").product(S.die(4), ["x", "y"], S.finite({ a: "1/6", b: "1/2", c: "1/3" })),
  "nested products": () => S.dice(2, 3).power(2),
};

describe("products are implicit but answer like their tables", () => {
  for (const [name, make] of Object.entries(cases)) {
    test(name, () => {
      const d = make();
      const t = listed(d);
      expect(d.factors).not.toBeNull();
      expect(t.factors).toBeNull();
      expect(d.count).toBe(t.count);
      expect(d.size).toBe(t.size);
      expect([...d.each()].map(([o, p]) => [JSON.stringify(o), str(p)])).toEqual(rows(t));
      expect(str(d)).toBe(str(t));
      const code = (o) => [...JSON.stringify(o)].reduce((s, c) => s + c.charCodeAt(0), 0);
      const f = (o) => code(o) % 3 === 0;
      const g = (o) => code(o) % 7;
      const h = (o) => (code(o) % 5) + JSON.stringify(o).charCodeAt(1);
      const some = new Set([t.outcomes[0], t.outcomes[4], t.outcomes[4], "nothing"]);
      for (const o of [t.outcomes[1], t.outcomes[t.size - 1], ["no"], "no"]) expect(str(d.p(o))).toBe(str(t.p(o)));
      expect(str(d.prob(f))).toBe(str(t.prob(f)));
      expect(str(d.prob(some))).toBe(str(t.prob(some)));
      const notFirst = (o) => JSON.stringify(o) !== JSON.stringify(t.outcomes[0]);
      expect(str(d.prob(f, notFirst))).toBe(str(t.prob(f, notFirst)));
      expect(str(d.prob(some, f))).toBe(str(t.prob(some, f)));
      expect(rows(d.given(f))).toEqual(rows(t.given(f)));
      expect(rows(d.given(some))).toEqual(rows(t.given(some)));
      expect(rows(d.given(t.outcomes[2]))).toEqual(rows(t.given(t.outcomes[2])));
      expect(rows(d.map(g))).toEqual(rows(t.map(g)));
      expect(rows(d.joint(g, h))).toEqual(rows(t.joint(g, h)));
      for (const which of [[0], [1], [0, 1], [1, 0], [0, 0]]) expect(rows(d.marginal(...which))).toEqual(rows(t.marginal(...which)));
      expect(d.independent(g, h)).toBe(t.independent(g, h));
      expect(str(d.mean(g))).toBe(str(t.mean(g)));
      expect(str(d.variance(g))).toBe(str(t.variance(g)));
      expect(str(d.moment(3, { of: g, central: true }))).toBe(str(t.moment(3, { of: g, central: true })));
      expect(str(d.covariance(g, h))).toBe(str(t.covariance(g, h)));
      expect(str(d.correlation(g, h, { bits: 24 }))).toBe(str(t.correlation(g, h, { bits: 24 })));
      expect(d.modes()).toEqual(t.modes());
      expect(d.samples(40, S.rng(7))).toEqual(t.samples(40, S.rng(7)));
      for (const ask of [(x) => x.mean(), (x) => x.cdf(2), (x) => x.sumOf(2), (x) => x.plus(S.die()), (x) => x.support, (x) => x.median()]) {
        let a, b;
        try {
          ask(d);
        } catch (e) {
          a = e.message;
        }
        try {
          ask(t);
        } catch (e) {
          b = e.message;
        }
        expect(a).toBeDefined();
        expect(a).toBe(b);
      }
      expect(() => d.given(() => false)).toThrow(/probability 0/);
      expect(() => d.prob(f, () => false)).toThrow(/probability 0/);
    });
  }
});

describe("large products", () => {
  test("nine dice stream without being listed", () => {
    const d = S.dice(9);
    expect(d.count).toBe(10077696n);
    expect(str(d.p([1, 2, 3, 4, 5, 6, 1, 2, 3]))).toBe("1/10077696");
    expect(str(d.p([1, 2, 3, 4, 5, 6, 1, 2, 7]))).toBe("0");
    // first equals last and some die is a six: 6⁸ − 5⁸ ... counted, then divided once
    expect(str(d.prob((r) => r[0] === r[8] && r.includes(6)))).toBe("1288991/10077696");
    expect(str(d.marginal(4))).toBe(str(S.die()));
    expect(d.marginal(0, 8).size).toBe(36);
    expect(d.sample(S.rng(1)).length).toBe(9);
    expect(d.toString()).toContain("10077696 outcomes");
    expect(() => d.entries).toThrow(/too many to list every outcome/);
    expect(() => d.table()).toThrow(/limit/);
    expect(() => d.modes()).toThrow(/most likely outcomes/);
  }, 30000);
  test("seven dice: Astra's numbers", () => {
    const d = S.dice(7);
    expect(str(d.prob((r) => r[0] === r[6] && r.includes(6)))).toBe("31031/279936");
    expect(str(d.sum().p(24))).toBe("24017/279936");
    expect(str(d.map((r) => r.reduce((a, b) => a + b, 0)).p(24))).toBe("24017/279936");
  });
  test("a sample reads the parts off one integer, as the table would", () => {
    const d = loaded().product(S.coin("1/3"), S.die(4));
    expect(d.samples(200, S.rng(3))).toEqual(listed(d).samples(200, S.rng(3)));
  });
  test("sum convolves the parts", () => {
    expect(str(S.dice(20).sum().p(70))).toBe(str(S.die().sumOf(20).p(70)));
    expect(str(S.dice(2).sum())).toBe(str(listed(S.dice(2)).sum()));
    expect(() => S.coin().product(S.die()).sum()).toThrow(/needs numbers/);
    expect(() => S.die().sum()).toThrow(/tuples/);
  });
  test("each() streams a table too", () => {
    expect([...S.coin().each()].length).toBe(2);
    expect(S.coin().count).toBe(2n);
  });
});
