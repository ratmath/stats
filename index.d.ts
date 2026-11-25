// Type declarations for @ratmath/stats.

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text ("1,200", or "600, plus 800 in graphics"). */
  estimate: string;
}

/** What the package is planned to contain. Nothing is built yet. */
export const plan: readonly PlannedPhase[];
