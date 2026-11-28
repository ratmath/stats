// Shared page furniture for the stats docs: the contents column, the theme
// switch and the code boxes. Add each topic page to PAGES as it is written.

import * as A from "./stats.js?v=4743369b49";
import { mountRunboxes } from "./runbox.js?v=4743369b49";

export { A };

export const PAGES = [
  { href: "index.html", title: "Overview" },
  { href: "reference.html", title: "Reference" },
  { href: "explore.html", title: "Explore" },
  { group: "Exact" },
  { href: "finite-probability.html", title: "Finite probability", blurb: "Sample spaces, conditioning, expectation and Bayes' rule in fractions" },
  { href: "many-dice.html", title: "Many dice", blurb: "Independent products kept implicit: stream, convolve and sample millions of rolls exactly" },
  { href: "discrete-distributions.html", title: "Discrete distributions", blurb: "Binomial, hypergeometric, geometric and Poisson probabilities" },
  { href: "descriptive.html", title: "Descriptive statistics", blurb: "Exact means, variances, named quantile rules, least squares" },
  { href: "exact-tests.html", title: "Exact tests", blurb: "Binomial, Fisher and permutation tests with p-values as fractions" },
  { group: "Certified" },
  { href: "continuous-distributions.html", title: "Continuous distributions", blurb: "Normal, t, chi-square, F, gamma, beta: enclosures to any precision" },
  { href: "densities.html", title: "Your own densities", blurb: "A distribution from a formula, functions of a random variable, order statistics" },
  { href: "certified-answers.html", title: "Certified answers", blurb: "Reading enclosures, and comparisons that can be too close to call" },
  { group: "Simulated" },
  { href: "sampling.html", title: "Seeded sampling", blurb: "A seeded generator, named streams, exact rational draws" },
  { href: "models.html", title: "Models with noise", blurb: "A known truth disturbed by noise, contamination, censoring and more" },
  { href: "experiments.html", title: "Experiments", blurb: "Coverage, size, power and bias, each with its sampling error" },
  { group: "Inference" },
  { href: "means.html", title: "Tests for means", blurb: "z and t tests, paired, Welch and pooled, with exact t² and certified p-values" },
  { href: "variances-proportions.html", title: "Variances and proportions", blurb: "Chi-square and F for spread; z tests and Wilson intervals for proportions" },
  { href: "chi-square.html", title: "Chi-square tests", blurb: "Goodness of fit and independence, with Yates' correction and Fisher's test beside them" },
  { href: "power.html", title: "Power and sample size", blurb: "Closed-form power of z tests, and the smallest n for a margin or a power" },
  { href: "estimation.html", title: "Fitting distributions", blurb: "Maximum likelihood and the method of moments for ten families" },
  { href: "testing-the-tests.html", title: "Testing the tests", blurb: "Seeded experiments check that intervals cover and tests keep their size" },
  { group: "Regression" },
  { href: "regression.html", title: "Linear regression", blurb: "Exact least squares, coefficient tests, prediction intervals, leverage, nested models" },
  { href: "anova.html", title: "Analysis of variance", blurb: "One-way ANOVA with exact sums of squares, as a regression, and two factors" },
  { group: "Nonparametric" },
  { href: "rank-tests.html", title: "Rank tests", blurb: "Sign, Wilcoxon, Mann–Whitney, Kruskal–Wallis, Spearman and Kendall, exact by counting" },
  { href: "kolmogorov-smirnov.html", title: "Kolmogorov–Smirnov", blurb: "One and two samples, with exact p-values by matrix powers and lattice paths" },
  { href: "bootstrap.html", title: "The bootstrap", blurb: "Exact over every resample for small samples, seeded otherwise" },
  { group: "Chains and beliefs" },
  { href: "markov-chains.html", title: "Markov chains", blurb: "Stationary distributions, classes, absorption and hitting times, exact" },
  { href: "bayesian-updating.html", title: "Bayesian updating", blurb: "Conjugate priors: exact posteriors and predictions, certified credible intervals" },
];

export const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};

// ---- chrome -------------------------------------------------------------------

function theme() {
  let saved = null;
  try {
    saved = localStorage.getItem("ratmath-theme");
  } catch {
    /* storage may be blocked */
  }
  const apply = (t) => {
    if (t === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    try {
      localStorage.setItem("ratmath-theme", t);
    } catch {
      /* ignore */
    }
    for (const b of document.querySelectorAll(".theme button")) b.classList.toggle("on", b.dataset.t === t);
  };
  if (saved) apply(saved);
  return h("div", { class: "theme", role: "group", "aria-label": "Theme" },
    ["auto", "light", "dark"].map((t) => h("button", { "data-t": t, class: (saved ?? "auto") === t ? "on" : "", onclick: () => apply(t) }, t)));
}

function chrome() {
  const here = location.pathname.split("/").pop() || "index.html";
  const aside = document.querySelector("aside");
  if (!aside) return;
  const nav = h("nav", { class: "toc", "aria-label": "Contents" });
  const sel = h("select", { id: "toc-select", "aria-label": "Page", onchange: (e) => (location.href = e.target.value) });
  let group = null;
  for (const p of PAGES) {
    if (p.group) {
      nav.append(h("div", { class: "group" }, p.group));
      group = h("optgroup", { label: p.group });
      sel.append(group);
      continue;
    }
    nav.append(h("a", { href: p.href, class: p.href === here ? "on" : null, "aria-current": p.href === here ? "page" : null }, p.title));
    (group ?? sel).append(h("option", { value: p.href, selected: p.href === here }, p.title));
  }
  aside.append(
    h("header", { class: "brand" }, h("a", { href: "index.html" }, h("div", { class: "mark" }, "RatMath ", h("span", {}, "stats"))), h("small", {}, "Probability and statistics")),
    sel, nav, theme(),
  );
}

/** Topic cards for the overview page (pages that have a blurb). */
export function topicCards(el) {
  el.append(...PAGES.filter((p) => p.blurb).map((p) => h("a", { href: p.href }, h("strong", {}, p.title), h("span", {}, p.blurb))));
}

chrome();
mountRunboxes(document, { scope: A });
