// Type declarations for @ratmath/stats.
import type { Rational, RationalInterval } from "@ratmath/core";
import type { Random } from "@ratmath/discrete";

export { rng, Random } from "@ratmath/discrete";

/** Anything read as an exact number: a JS number is read by its decimal form. */
export type Num = Rational | number | bigint | string;
/** An exact value or a certified enclosure of one. */
export type Value = Rational | RationalInterval;
/** A datum: a number, or an interval as [low, high], "1.2:1.4" or "1.5[+-0.1]". */
export type Datum = Num | RationalInterval | [Num, Num];
export type Seed = number | string | bigint;
export interface Bits { bits?: number }
export type Alternative = "two-sided" | "greater" | "less";

export class StatsError extends Error {
  code: string;
}
export function toQ(v: Num | RationalInterval, what?: string): Rational;
export function toValue(v: Datum, what?: string): Value;
/** Decimal text; for an enclosure, the digits both ends share, then "…". */
export function show(v: Value | number, digits?: number): string;
export function compare(value: Value, level: Num): "below" | "above" | "equal" | "too close to call";
export function sqrt(v: Datum, bits?: number): Value;
export function toNumber(v: Value | number): number;

// ---- finite probability -------------------------------------------------------

export type Event<T> = ((outcome: T) => boolean) | Set<T> | T;
export class Finite<T = unknown> {
  constructor(entries: Iterable<[T, Num]>, options?: { normalize?: boolean });
  readonly entries: readonly (readonly [T, Rational])[];
  readonly size: number;
  readonly outcomes: T[];
  readonly support: T[];
  p(outcome: T): Rational;
  prob(event: Event<T>, given?: Event<T>): Rational;
  given(event: Event<T>): Finite<T>;
  filter(event: Event<T>): Finite<T>;
  map<U>(f: (outcome: T) => U): Finite<U>;
  product(...others: Finite<any>[]): Finite<any[]>;
  power(n: number): Finite<T[]>;
  plus(other: Finite<any>): Finite<Rational>;
  sumOf(n: number): Finite<Rational>;
  marginal(...which: number[]): Finite<any>;
  joint(...fs: ((outcome: T) => unknown)[]): Finite<unknown[]>;
  independent(a: Event<T>, b: Event<T>): boolean;
  mean(f?: (outcome: T) => Num): Rational;
  expect(f: (outcome: T) => Num): Rational;
  moment(k: number, options?: { central?: boolean; of?: (outcome: T) => Num }): Rational;
  variance(f?: (outcome: T) => Num): Rational;
  sd(f?: (outcome: T) => Num, options?: Bits): Value;
  covariance(f: (outcome: T) => Num, g: (outcome: T) => Num): Rational;
  correlation(f: (outcome: T) => Num, g: (outcome: T) => Num, options?: Bits): Value;
  cdf(x: Num): Rational;
  quantile(p: Num): T;
  median(): T;
  modes(): T[];
  sample(rng: Random): T;
  samples(n: number, rng: Random): T[];
  table(options?: { sort?: "probability" | "outcome" }): { outcome: T; probability: Rational }[];
  toString(): string;
}
export function finite<T>(outcomes: T[], weights?: Num[], options?: { normalize?: boolean }): Finite<T>;
export function finite<T>(spec: Map<T, Num> | Iterable<T>, options?: { normalize?: boolean }): Finite<T>;
export function finite(spec: Record<string, Num>, options?: { normalize?: boolean }): Finite<string>;
export function die(sides?: number): Finite<number>;
export function dice(n: number, sides?: number): Finite<number[]>;
export function coin(p?: Num): Finite<"H" | "T">;
/** The 52 cards as text such as "Q♥". */
export function deck(): string[];
export function draw<T>(items: T[], k: number, options?: { replace?: boolean; ordered?: boolean }): Finite<T[]>;
export function bayes<T>(
  prior: Finite<T> | Record<string, Num> | Map<T, Num>,
  likelihood: ((hypothesis: T) => Num) | Record<string, Num> | Map<T, Num>,
): { posterior: Finite<T>; evidence: Rational; table: { hypothesis: T; prior: Rational; likelihood: Rational; joint: Rational; posterior: Rational }[] };

// ---- discrete distributions ---------------------------------------------------

export class Discrete {
  readonly kind: "exact" | "certified";
  readonly name: string;
  readonly min: number;
  /** null when the values are unbounded. */
  readonly max: number | null;
  readonly mean: Rational;
  readonly variance: Rational;
  readonly mode: number;
  pmf(k: Num, options?: Bits): Value;
  cdf(k: Num, options?: Bits): Value;
  sf(k: Num, options?: Bits): Value;
  /** P(a ≤ X ≤ b). */
  prob(a: Num, b: Num, options?: Bits): Value;
  /** The least k with P(X ≤ k) ≥ p. */
  quantile(p: Num, options?: Bits): number;
  median(options?: Bits): number;
  sd(options?: Bits): Value;
  sample(rng: Random): number;
  samples(n: number, rng: Random): number[];
  toFinite(): Finite<number>;
  toString(): string;
}
export function bernoulli(p: Num): Discrete;
export function binomial(n: number, p: Num): Discrete;
export function geometric(p: Num, options?: { failures?: boolean }): Discrete;
export function negativeBinomial(r: number, p: Num): Discrete;
export function hypergeometric(N: number, K: number, n: number): Discrete;
export function discreteUniform(a: number, b: number): Discrete;
export function poisson(lambda: Num): Discrete;

// ---- descriptive statistics ---------------------------------------------------

export type QuantileMethod = "linear" | "weibull" | "hazen" | "median-unbiased" | "inverse" | "lower" | "higher" | "nearest" | "midpoint";
export function count(data: Datum[]): number;
export function sum(data: Datum[]): Value;
export function mean(data: Datum[]): Value;
export function weightedMean(data: Datum[], weights: Num[]): Value;
export function variance(data: Datum[], options?: { population?: boolean }): Value;
export function sd(data: Datum[], options?: { population?: boolean; bits?: number }): Value;
export function moment(data: Datum[], k: number, options?: { central?: boolean }): Value;
export function skewness(data: Num[], options?: Bits): Value;
export function kurtosis(data: Num[], options?: { excess?: boolean }): Rational;
export function min(data: Datum[]): Value;
export function max(data: Datum[]): Value;
export function range(data: Datum[]): Value;
export function quantile(data: Num[], p: Num, options?: { method?: QuantileMethod }): Rational;
export function median(data: Num[]): Rational;
/** [q1, median, q3]. */
export function quartiles(data: Num[], options?: { method?: QuantileMethod }): [Rational, Rational, Rational];
export function iqr(data: Num[], options?: { method?: QuantileMethod }): Rational;
/** [min, q1, median, q3, max]. */
export function fiveNumber(data: Num[], options?: { method?: QuantileMethod }): [Rational, Rational, Rational, Rational, Rational];
export function mad(data: Num[]): Rational;
export function trimmedMean(data: Num[], fraction: Num): Rational;
export function frequencies(data: Num[]): { value: Rational; count: number; proportion: Rational; cumulative: Rational }[];
export function modes(data: Num[]): Rational[];
export function histogram(
  data: Num[],
  options?: { edges?: Num[]; width?: Num; start?: Num; bins?: number },
): { from: Rational; to: Rational; count: number; proportion: Rational; density: Rational }[];
export function covariance(x: Datum[], y: Datum[], options?: { population?: boolean }): Value;
export function correlation(x: Num[], y: Num[], options?: Bits): Value;
export interface Fit {
  slope: Rational;
  intercept: Rational;
  r2: Rational;
  residuals: Rational[];
  fitted: Rational[];
  sse: Rational;
  sst: Rational;
  n: number;
  predict(x: Num): Rational;
}
export function leastSquares(x: Num[], y: Num[]): Fit;
export function summary(data: Datum[], options?: Bits): {
  n: number; mean: Value; min: Value; q1: Value; median: Value; q3: Value; max: Value; iqr: Value; range: Value;
  variance?: Value; sd?: Value; outliers?: Rational[]; modes?: Rational[];
};
/** The same statistics in floating point on plain numbers, for simulations. */
export const fast: {
  sum(data: number[]): number;
  mean(data: number[]): number;
  variance(data: number[], options?: { population?: boolean }): number;
  sd(data: number[], options?: { population?: boolean }): number;
  quantile(data: number[], p: number): number;
  median(data: number[]): number;
  min(data: number[]): number;
  max(data: number[]): number;
  trimmedMean(data: number[], fraction: number): number;
  covariance(x: number[], y: number[]): number;
  correlation(x: number[], y: number[]): number;
  leastSquares(x: number[], y: number[]): { slope: number; intercept: number; r2: number; n: number };
};

// ---- continuous distributions -------------------------------------------------

export class Continuous {
  readonly kind: "certified";
  readonly name: string;
  readonly params: Record<string, Rational>;
  /** The ends of the support; null for infinite. */
  readonly low: Rational | null;
  readonly high: Rational | null;
  /** null where the moment does not exist. */
  readonly mean: Value | null;
  readonly variance: Value | null;
  pdf(x: Num | RationalInterval, options?: Bits): Value;
  cdf(x: Num | RationalInterval, options?: Bits): Value;
  sf(x: Num | RationalInterval, options?: Bits): Value;
  prob(a: Num, b: Num, options?: Bits): Value;
  quantile(p: Num, options?: Bits): Value;
  median(options?: Bits): Value;
  sd(options?: Bits): Value | null;
  sample(rng: Random): number;
  samples(n: number, rng: Random): number[];
  sampleExact(rng: Random, options?: Bits): Rational;
  samplesExact(n: number, rng: Random, options?: Bits): Rational[];
  readonly fast: { pdf(x: number): number; cdf(x: number): number; quantile(p: number): number };
  toString(): string;
}
export function normal(mu?: Num, sigma?: Num): Continuous;
export function lognormal(mu?: Num, sigma?: Num): Continuous;
export function gamma(shape: Num, rate?: Num): Continuous;
export function chiSquare(k: Num): Continuous;
export function exponential(rate?: Num): Continuous;
export function weibull(shape: Num, scale?: Num): Continuous;
export function beta(a: Num, b: Num): Continuous;
export function studentT(df: Num): Continuous;
export function fDist(d1: Num, d2: Num): Continuous;
export function uniform(a?: Num, b?: Num): Continuous;
export function cauchy(location?: Num, scale?: Num): Continuous;

// ---- exact tests ----------------------------------------------------------------

export class TestResult {
  readonly test: string;
  readonly statistic: unknown;
  readonly pValue: Value;
  readonly alternative: Alternative;
  readonly kind: "exact" | "certified" | "simulated";
  decision(alpha?: Num): "reject" | "keep" | "too close to call";
  /** Refuses when the p-value cannot be told apart from alpha. */
  reject(alpha?: Num): boolean;
  toString(): string;
  [field: string]: unknown;
}
export interface ProportionInterval { low: Value; high: Value; level: Rational; method: "exact" | "wilson" | "wald"; estimate: Rational }
export function proportionInterval(successes: number, trials: number, options?: { level?: Num; method?: "exact" | "wilson" | "wald"; bits?: number }): ProportionInterval;
export function binomialTest(successes: number, trials: number, p0?: Num, options?: { alternative?: Alternative }): TestResult & { interval(options?: { level?: Num; method?: "exact" | "wilson" | "wald"; bits?: number }): ProportionInterval };
export function fisherTest(table: [[number, number], [number, number]], options?: { alternative?: Alternative }): TestResult & { oddsRatio: Rational | number | null };
export function permutationTest(
  x: Num[],
  y: Num[],
  options?: { statistic?: "mean" | "median" | ((values: Rational[]) => Num); alternative?: Alternative; paired?: boolean; limit?: number; resamples?: number; seed?: Seed },
): TestResult & { arrangements: bigint; pValueInterval?: { low: Value; high: Value; level: Rational }; seed?: Seed };

// ---- models and experiments -------------------------------------------------------

export type Source = Continuous | Discrete | Finite<any> | ((rng: Random) => Num) | Num;
export interface ModelSpec {
  /** A number, the coefficients [b0, b1, …] of a polynomial in x, or a function of x. */
  truth: Num | Num[] | ((x: any) => Num);
  /** A list reused in order, a distribution, or a function of the generator. */
  x?: Num[] | Source;
  noise?: Source;
  /** Multiplies the noise: unequal spread. */
  spread?: (x: any) => Num;
  contamination?: { rate: Num; noise?: Source; value?: Source };
  /** Added to the recorded x; y is made from the true x. */
  xNoise?: Source;
  round?: Num | { x?: Num; y?: Num };
  censor?: { below?: Num; above?: Num };
  /** The probability that an observation is lost. */
  missing?: Num;
  name?: string;
}
export class Dataset<N = number> {
  readonly kind: "simulated";
  readonly x: N[] | undefined;
  readonly y: N[];
  readonly trueX: N[] | undefined;
  readonly trueY: N[];
  readonly contaminated: boolean[];
  readonly censored: boolean[];
  readonly missing: boolean[];
  readonly seed: Seed;
  readonly exact: boolean;
  readonly model: string;
  readonly n: number;
  complete(): { x: N[] | undefined; y: N[] };
  toString(): string;
}
export class Model {
  constructor(spec: ModelSpec);
  readonly spec: Readonly<ModelSpec>;
  with(changes: Partial<ModelSpec>): Model;
  truthAt(x?: Num, exact?: boolean): number | Rational;
  sample(n: number, seed?: Seed | Random, options?: { exact?: false; bits?: number }): Dataset<number>;
  sample(n: number, seed: Seed | Random, options: { exact: true; bits?: number }): Dataset<Rational | null>;
  toString(): string;
}
export function model(spec: ModelSpec): Model;

export type DataSource<D> = { model: Model; n: number } | { draw: (rng: Random, run: number) => D };
export interface RunOptions { runs?: number; seed?: Seed }
export class SimulatedRate {
  readonly kind: "simulated";
  readonly rate: Rational;
  readonly hits: number;
  readonly runs: number;
  readonly seed: Seed;
  readonly interval: { low: Value; high: Value; level: Rational };
  readonly truth?: Rational;
  readonly meanWidth?: number;
  readonly alpha?: Rational;
  readonly effect?: unknown;
  /** Where a stated value lies against the interval. */
  compare(value: Num): "below" | "inside" | "above";
  toString(): string;
}
export function howOften<D = Dataset>(options: DataSource<D> & RunOptions & { event: (data: D, run: number) => boolean; level?: Num }): SimulatedRate;
export function coverage<D = Dataset>(
  options: DataSource<D> & RunOptions & { interval: (data: D) => [Datum, Datum] | { low: Datum; high: Datum }; truth: Num; level?: Num },
): SimulatedRate;
export function rejectionRate<D = Dataset>(options: DataSource<D> & RunOptions & { test: (data: D) => number | Value | TestResult | boolean; alpha?: Num; level?: Num }): SimulatedRate;
export function powerCurve<E, D = Dataset>(
  options: RunOptions & { effects: E[]; model?: (effect: E) => Model; n?: number; draw?: (effect: E) => (rng: Random, run: number) => D; test: (data: D) => number | Value | TestResult | boolean; alpha?: Num; level?: Num },
): SimulatedRate[];
export class SimulatedEstimates {
  readonly kind: "simulated";
  readonly estimates: number[];
  readonly mean: number;
  readonly sd: number;
  /** Of the mean of the estimates. */
  readonly standardError: number;
  readonly runs: number;
  readonly seed: Seed;
  readonly truth?: number;
  readonly bias?: number;
  readonly rmse?: number;
  toString(): string;
}
export function estimatorStudy<D = Dataset>(options: DataSource<D> & RunOptions & { estimator: (data: D, run: number) => number | Value; truth?: Num }): SimulatedEstimates;

// ---- what is still planned --------------------------------------------------------

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text ("1,200", or "800 in graphics"). */
  estimate: string;
}
/** The phases that are not built yet. */
export const plan: readonly PlannedPhase[];
